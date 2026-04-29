#include "bcl_parser_cuda.h"

#include <cuda_runtime.h>

#include <cstdio>
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

__global__ void decode_bcl_kernel(const char* const* d_bcl_data,
                                  const int* d_read_structure,
                                  char* d_seq,
                                  char* d_qual,
                                  int num_cycles,
                                  int total_seq_len,
                                  size_t batch_size,
                                  int r1_len,
                                  int i1_len,
                                  int i2_len) {
    size_t idx = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= batch_size) return;

    static const char bases[] = {'A', 'C', 'G', 'T'};

    int pos_r1 = 0, pos_i1 = 0, pos_i2 = 0, pos_r2 = 0;
    size_t cluster_offset = idx * static_cast<size_t>(total_seq_len);

    for (int cycle = 0; cycle < num_cycles; ++cycle) {
        unsigned char bcl_byte = static_cast<unsigned char>(d_bcl_data[cycle][idx]);
        int base_code = bcl_byte & 0x03;
        int quality_val = (bcl_byte >> 2) & 0x3F;
        char base = bases[base_code];
        char quality_char = static_cast<char>(quality_val + 33);

        int seg = d_read_structure[cycle];
        size_t out;
        if (seg == 0) {
            out = cluster_offset + pos_r1++;
        } else if (seg == 1) {
            out = cluster_offset + r1_len + pos_i1++;
        } else if (seg == 2) {
            out = cluster_offset + r1_len + i1_len + pos_i2++;
        } else if (seg == 3) {
            out = cluster_offset + r1_len + i1_len + i2_len + pos_r2++;
        } else {
            continue;
        }
        d_seq[out] = base;
        d_qual[out] = quality_char;
    }
}

}  // namespace

struct CudaDecodeContext {
    cudaStream_t stream = nullptr;
    int num_cycles = 0;
    int* d_read_structure = nullptr;
    const char** d_bcl_ptrs_dev = nullptr;
    std::vector<char*> d_cycle_buffers;        // device pointers, one per cycle
    std::vector<const char*> h_cycle_ptrs;     // host-side mirror for the above
    size_t per_cycle_capacity = 0;
};

CudaDecodeContext* decode_context_create(const LaneBclData& lane) {
    auto* ctx = new CudaDecodeContext();
    try {
        CUDA_CHECK(cudaStreamCreate(&ctx->stream));
        ctx->num_cycles = lane.total_cycles;
        CUDA_CHECK(cudaMalloc(&ctx->d_read_structure, ctx->num_cycles * sizeof(int)));
        CUDA_CHECK(cudaMemcpyAsync(ctx->d_read_structure,
                                   lane.read_segments.data(),
                                   ctx->num_cycles * sizeof(int),
                                   cudaMemcpyHostToDevice,
                                   ctx->stream));
        CUDA_CHECK(cudaMalloc(&ctx->d_bcl_ptrs_dev, ctx->num_cycles * sizeof(char*)));
        ctx->d_cycle_buffers.assign(ctx->num_cycles, nullptr);
        ctx->h_cycle_ptrs.assign(ctx->num_cycles, nullptr);
        CUDA_CHECK(cudaStreamSynchronize(ctx->stream));
    } catch (...) {
        decode_context_destroy(ctx);
        throw;
    }
    return ctx;
}

void decode_context_destroy(CudaDecodeContext* ctx) {
    if (!ctx) return;
    for (auto* p : ctx->d_cycle_buffers) {
        if (p) cudaFree(p);
    }
    if (ctx->d_bcl_ptrs_dev) cudaFree(ctx->d_bcl_ptrs_dev);
    if (ctx->d_read_structure) cudaFree(ctx->d_read_structure);
    if (ctx->stream) cudaStreamDestroy(ctx->stream);
    delete ctx;
}

bool decode_bcl_batch(CudaDecodeContext* ctx,
                      const LaneBclData& lane,
                      size_t batch_start,
                      size_t batch_size,
                      char* d_seq,
                      char* d_qual) {
    if (batch_size == 0) return true;
    if (batch_start + batch_size > lane.num_clusters) {
        throw std::runtime_error("decode_bcl_batch: batch out of range");
    }

    if (ctx->per_cycle_capacity < batch_size) {
        for (auto*& p : ctx->d_cycle_buffers) {
            if (p) {
                cudaFree(p);
                p = nullptr;
            }
        }
        ctx->per_cycle_capacity = 0;
        for (int c = 0; c < ctx->num_cycles; ++c) {
            cudaError_t e = cudaMalloc(&ctx->d_cycle_buffers[c], batch_size);
            if (e != cudaSuccess) {
                for (auto*& p : ctx->d_cycle_buffers) {
                    if (p) {
                        cudaFree(p);
                        p = nullptr;
                    }
                }
                return false;
            }
            ctx->h_cycle_ptrs[c] = ctx->d_cycle_buffers[c];
        }
        ctx->per_cycle_capacity = batch_size;
        CUDA_CHECK(cudaMemcpyAsync(ctx->d_bcl_ptrs_dev,
                                   ctx->h_cycle_ptrs.data(),
                                   ctx->num_cycles * sizeof(char*),
                                   cudaMemcpyHostToDevice,
                                   ctx->stream));
    }

    const int total_seq_len = lane.r1_len + lane.i1_len + lane.i2_len + lane.r2_len;
    if (total_seq_len <= 0) {
        return true;
    }

    CUDA_CHECK(cudaMemsetAsync(d_seq, 'N', batch_size * total_seq_len, ctx->stream));
    CUDA_CHECK(cudaMemsetAsync(d_qual, '!', batch_size * total_seq_len, ctx->stream));

    for (int c = 0; c < ctx->num_cycles; ++c) {
        const uint8_t* src = lane.bcl[c].data() + batch_start;
        CUDA_CHECK(cudaMemcpyAsync(ctx->d_cycle_buffers[c],
                                   src,
                                   batch_size,
                                   cudaMemcpyHostToDevice,
                                   ctx->stream));
    }

    constexpr int kThreads = 256;
    int blocks = static_cast<int>((batch_size + kThreads - 1) / kThreads);
    decode_bcl_kernel<<<blocks, kThreads, 0, ctx->stream>>>(
        ctx->d_bcl_ptrs_dev,
        ctx->d_read_structure,
        d_seq,
        d_qual,
        ctx->num_cycles,
        total_seq_len,
        batch_size,
        lane.r1_len,
        lane.i1_len,
        lane.i2_len);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(ctx->stream));
    return true;
}
