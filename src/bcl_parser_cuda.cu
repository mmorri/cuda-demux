#include "bcl_parser_cuda.h"

#include <cuda_runtime.h>
#include <omp.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

#define CUDA_CHECK(expr)                                                                       \
    do {                                                                                       \
        cudaError_t _e = (expr);                                                               \
        if (_e != cudaSuccess) {                                                               \
            throw std::runtime_error(std::string("CUDA error at " __FILE__ ":") +              \
                                     std::to_string(__LINE__) + ": " + cudaGetErrorString(_e));\
        }                                                                                      \
    } while (0)

constexpr int kTile = 32;
constexpr int kTileRows = 8;

// Transposes cycle-major BCL bytes into cluster-major sequence/quality text.
//
// Input  d_cycles: [num_cycles][batch_stride] bytes, (q << 2) | base, 0 = no-call.
// Output d_seq/d_qual: [batch][total_seq_len] chars; cycle c lands in column
// d_out_col[c] (its position after regrouping cycles into R1, I1, I2, R2).
//
// A block handles a 32-cycle x 32-cluster tile through shared memory so both the
// global reads (along clusters) and the global writes (along a cluster's row)
// are contiguous, instead of every thread striding by total_seq_len.
__global__ void decode_bcl_kernel(const unsigned char* __restrict__ d_cycles,
                                  size_t batch_stride,
                                  const int* __restrict__ d_out_col,
                                  char* __restrict__ d_seq,
                                  char* __restrict__ d_qual,
                                  int num_cycles,
                                  int total_seq_len,
                                  size_t batch_size) {
    __shared__ unsigned char tile[kTile][kTile + 1];

    const int cyc0 = blockIdx.y * kTile;
    const size_t clu0 = static_cast<size_t>(blockIdx.x) * kTile;

    // Load: threadIdx.x walks clusters (contiguous in memory for a cycle).
    for (int cy = threadIdx.y; cy < kTile; cy += kTileRows) {
        const int c = cyc0 + cy;
        const size_t k = clu0 + threadIdx.x;
        tile[cy][threadIdx.x] = (c < num_cycles && k < batch_size)
                                    ? d_cycles[static_cast<size_t>(c) * batch_stride + k]
                                    : 0;
    }
    __syncthreads();

    // Store: threadIdx.x walks cycles (contiguous columns in a cluster's row).
    const int c = cyc0 + threadIdx.x;
    if (c >= num_cycles) return;
    const int col = d_out_col[c];
    for (int cl = threadIdx.y; cl < kTile; cl += kTileRows) {
        const size_t k = clu0 + cl;
        if (k >= batch_size) break;
        const unsigned char b = tile[threadIdx.x][cl];
        const int q = b >> 2;
        const size_t out = k * static_cast<size_t>(total_seq_len) + col;
        // Q0 is a no-call: emit 'N' with Q2 ('#'), matching bcl2fastq/bcl-convert.
        d_seq[out] = q ? "ACGT"[b & 3] : 'N';
        d_qual[out] = q ? static_cast<char>(q + 33) : '#';
    }
}

struct Slot {
    unsigned char* h_stage = nullptr;   // pinned, [num_cycles][capacity]
    unsigned char* d_cycles = nullptr;  // device, same layout
};

}  // namespace

struct CudaDecodeContext {
    cudaStream_t stream = nullptr;
    int num_cycles = 0;
    int total_seq_len = 0;
    size_t capacity = 0;
    int* d_out_col = nullptr;
    std::vector<Slot> slots;
};

cudaStream_t decode_context_stream(CudaDecodeContext* ctx) { return ctx->stream; }

CudaDecodeContext* decode_context_create(int num_cycles, size_t batch_capacity, int num_slots) {
    auto* ctx = new CudaDecodeContext();
    try {
        CUDA_CHECK(cudaStreamCreate(&ctx->stream));
        ctx->num_cycles = num_cycles;
        ctx->capacity = batch_capacity;
        CUDA_CHECK(cudaMalloc(&ctx->d_out_col, num_cycles * sizeof(int)));
        ctx->slots.resize(num_slots);
        const size_t bytes = static_cast<size_t>(num_cycles) * batch_capacity;
        for (Slot& s : ctx->slots) {
            CUDA_CHECK(cudaHostAlloc((void**)&s.h_stage, bytes, cudaHostAllocDefault));
            CUDA_CHECK(cudaMalloc(&s.d_cycles, bytes));
        }
    } catch (...) {
        decode_context_destroy(ctx);
        throw;
    }
    return ctx;
}

void decode_context_set_lane(CudaDecodeContext* ctx, const LaneBclData& lane) {
    if (lane.total_cycles != ctx->num_cycles ||
        static_cast<int>(lane.read_segments.size()) != lane.total_cycles) {
        throw std::runtime_error("decode_context_set_lane: lane cycle structure does not match context");
    }
    ctx->total_seq_len = lane.r1_len + lane.i1_len + lane.i2_len + lane.r2_len;

    // Output column of each cycle: cycles of a segment are contiguous in the
    // segment's slice of the row, and segments are laid out R1 | I1 | I2 | R2.
    const int seg_base[4] = {0, lane.r1_len, lane.r1_len + lane.i1_len,
                             lane.r1_len + lane.i1_len + lane.i2_len};
    int seg_pos[4] = {0, 0, 0, 0};
    std::vector<int> out_col(ctx->num_cycles, 0);
    for (int c = 0; c < ctx->num_cycles; ++c) {
        const int seg = lane.read_segments[c];
        if (seg < 0 || seg > 3) throw std::runtime_error("invalid read segment");
        out_col[c] = seg_base[seg] + seg_pos[seg]++;
    }
    CUDA_CHECK(cudaStreamSynchronize(ctx->stream));
    CUDA_CHECK(cudaMemcpy(ctx->d_out_col, out_col.data(), ctx->num_cycles * sizeof(int),
                          cudaMemcpyHostToDevice));
}

void decode_context_destroy(CudaDecodeContext* ctx) {
    if (!ctx) return;
    if (ctx->stream) cudaStreamSynchronize(ctx->stream);
    for (Slot& s : ctx->slots) {
        if (s.h_stage) cudaFreeHost(s.h_stage);
        if (s.d_cycles) cudaFree(s.d_cycles);
    }
    if (ctx->d_out_col) cudaFree(ctx->d_out_col);
    if (ctx->stream) cudaStreamDestroy(ctx->stream);
    delete ctx;
}

bool decode_bcl_batch(CudaDecodeContext* ctx,
                      const LaneBclData& lane,
                      size_t batch_start,
                      size_t batch_size,
                      int slot_idx,
                      char* d_seq,
                      char* d_qual) {
    if (batch_size == 0 || ctx->total_seq_len <= 0) return true;
    if (lane.total_cycles != ctx->num_cycles) {
        throw std::runtime_error("decode_bcl_batch: call decode_context_set_lane first");
    }
    if (batch_start + batch_size > lane.num_clusters) {
        throw std::runtime_error("decode_bcl_batch: batch out of range");
    }
    if (batch_size > ctx->capacity) {
        throw std::runtime_error("decode_bcl_batch: batch exceeds context capacity");
    }
    if (slot_idx < 0 || slot_idx >= static_cast<int>(ctx->slots.size())) {
        throw std::runtime_error("decode_bcl_batch: bad slot");
    }
    Slot& slot = ctx->slots[slot_idx];
    const int num_cycles = ctx->num_cycles;

    // Gather this batch's slice of every cycle into pinned memory so the upload
    // is one asynchronous DMA instead of hundreds of pageable copies.
    #pragma omp parallel for schedule(static)
    for (int c = 0; c < num_cycles; ++c) {
        std::memcpy(slot.h_stage + static_cast<size_t>(c) * batch_size,
                    lane.bcl[c].data() + batch_start, batch_size);
    }
    CUDA_CHECK(cudaMemcpyAsync(slot.d_cycles, slot.h_stage,
                               static_cast<size_t>(num_cycles) * batch_size,
                               cudaMemcpyHostToDevice, ctx->stream));

    const dim3 block(kTile, kTileRows);
    const dim3 grid(static_cast<unsigned>((batch_size + kTile - 1) / kTile),
                    static_cast<unsigned>((num_cycles + kTile - 1) / kTile));
    decode_bcl_kernel<<<grid, block, 0, ctx->stream>>>(
        slot.d_cycles, batch_size, ctx->d_out_col, d_seq, d_qual,
        num_cycles, ctx->total_seq_len, batch_size);
    CUDA_CHECK(cudaGetLastError());
    return true;
}
