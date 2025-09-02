#include "bcl_parser_cuda.h"
#include <cuda_runtime.h>
#include <iostream>
#include <numeric>
#include <cstdlib>

// Helper macro for CUDA error checking
#define CUDA_CHECK(err) { \
    cudaError_t err_ = (err); \
    if (err_ != cudaSuccess) { \
        std::cerr << "CUDA error in " << __FILE__ << " at line " << __LINE__ << ": " << cudaGetErrorString(err_) << std::endl; \
        exit(EXIT_FAILURE); \
    } \
}

// CUDA kernel to decode BCL data for all clusters in parallel
__global__ void decode_bcl_kernel(
    const char** d_bcl_data,      // Device pointers to raw BCL data for each cycle
    const int* d_read_structure,  // Defines which cycle belongs to which read segment
    char* d_output_sequences,     // Output buffer for all decoded sequences concatenated
    char* d_output_qualities,     // Output buffer for all decoded quality scores concatenated
    int num_cycles,               // Total number of cycles
    int total_sequence_length,    // Total length of all read segments combined
    size_t num_clusters,          // Total number of clusters (reads)
    int read1_len,                // R1 cycles count
    int index1_len,               // I1 cycles count
    int index2_len,               // I2 cycles count
    int read2_len                 // R2 cycles count
) {
    // Use 64-bit indexing to avoid overflow for large batches
    size_t cluster_idx = static_cast<size_t>(blockIdx.x) * static_cast<size_t>(blockDim.x) + static_cast<size_t>(threadIdx.x);
    if (cluster_idx >= num_clusters) {
        return;
    }

    static const char bases[] = {'A', 'C', 'G', 'T', 'N'};

    // Per-segment write positions to enforce canonical layout: R1, I1, I2, R2
    int pos_r1 = 0, pos_i1 = 0, pos_i2 = 0, pos_r2 = 0;
    size_t base_idx = cluster_idx * static_cast<size_t>(total_sequence_length);
    for (int cycle = 0; cycle < num_cycles; ++cycle) {
        // Each BCL file buffer contains data for this batch of clusters for one cycle
        // cluster_idx is local to this batch (0 to num_clusters-1)
        char bcl_byte = d_bcl_data[cycle][cluster_idx];

        // Decode base and quality
        // BCL encoding: bits 0-1: base (0=A, 1=C, 2=G, 3=T), bits 2-7: quality
        int base_idx = bcl_byte & 0x03;
        int quality_val = (bcl_byte >> 2) & 0x3F;

        char base = (quality_val == 0) ? 'N' : bases[base_idx];
        char quality_char = static_cast<char>(quality_val + 33); // Phred+33

        // Determine where to write the decoded base and quality
        int read_segment = d_read_structure[cycle];
        if (read_segment >= 0) { // Negative values can be used to skip cycles if needed
            size_t idx = base_idx; // start of this cluster's concatenated output
            if (read_segment == 0) { // R1
                idx += static_cast<size_t>(pos_r1++);
            } else if (read_segment == 1) { // I1
                idx += static_cast<size_t>(read1_len + pos_i1++);
            } else if (read_segment == 2) { // I2
                idx += static_cast<size_t>(read1_len + index1_len + pos_i2++);
            } else if (read_segment == 3) { // R2
                idx += static_cast<size_t>(read1_len + index1_len + index2_len + pos_r2++);
            } else {
                // Unknown segment; skip
                continue;
            }
            d_output_sequences[idx] = base;
            d_output_qualities[idx] = quality_char;
        }
    }
}

// Host function to manage BCL decoding on the GPU
void decode_bcl_data_cuda(
    const std::vector<char*>& h_bcl_data,
    const std::vector<size_t>& h_bcl_sizes,
    const std::vector<int>& h_read_structure,
    std::vector<Read>& reads,
    size_t num_clusters
) {
    if (h_bcl_data.empty() || num_clusters == 0) {
        return;
    }

    int num_cycles = h_bcl_data.size();
    
    // Calculate actual lengths of each segment first
    int read1_len = 0, index1_len = 0, index2_len = 0, read2_len = 0;
    for (int seg : h_read_structure) {
        if (seg == 0) read1_len++;
        else if (seg == 1) index1_len++;
        else if (seg == 2) index2_len++;
        else if (seg == 3) read2_len++;
    }
    
    // Total sequence length is the sum of all segments
    int total_sequence_length = read1_len + index1_len + index2_len + read2_len;
    
    std::cout << "Debug: num_cycles=" << num_cycles 
              << ", total_sequence_length=" << total_sequence_length 
              << " (R1=" << read1_len << ", I1=" << index1_len 
              << ", I2=" << index2_len << ", R2=" << read2_len << ")" << std::endl;
    std::cout << "Debug: h_read_structure size=" << h_read_structure.size() << std::endl;

    // --- Streaming approach: process clusters in GPU-sized batches ---
    // Calculate a safe batch size based on available GPU memory
    size_t free_mem = 0, total_mem = 0;
    // Prefer cudaMemGetInfo; if unavailable, fall back to device properties (totalGlobalMem)
    cudaError_t meminfo_err = cudaMemGetInfo(&free_mem, &total_mem);
    if (meminfo_err != cudaSuccess) {
        int dev = 0;
        cudaDeviceProp prop{};
        cudaError_t dev_err = cudaGetDevice(&dev);
        cudaError_t prop_err = (dev_err == cudaSuccess) ? cudaGetDeviceProperties(&prop, dev) : cudaErrorUnknown;
        if (dev_err == cudaSuccess && prop_err == cudaSuccess && prop.totalGlobalMem > 0) {
            // Allow override of usable fraction via env var
            double frac = 0.60; // default to 60% of total mem for working set
            if (const char* env = std::getenv("CUDA_DEMUX_MEM_FRACTION")) {
                try {
                    double v = std::stod(env);
                    if (v > 0.05 && v <= 0.95) frac = v;
                } catch (...) {}
            }
            total_mem = static_cast<size_t>(prop.totalGlobalMem);
            free_mem  = static_cast<size_t>(prop.totalGlobalMem * frac);
            std::cerr << "Info: cudaMemGetInfo unavailable (" << cudaGetErrorString(meminfo_err)
                      << "); using totalGlobalMem with fraction " << frac << "." << std::endl;
        } else {
            // Final fallback: modest static budget
            total_mem = (size_t)3ULL * 1024ULL * 1024ULL * 1024ULL;
            free_mem  = (size_t)2ULL * 1024ULL * 1024ULL * 1024ULL;
            std::cerr << "Warning: Could not query device properties; using static memory estimates." << std::endl;
        }
    }

    // Bytes per cluster on device: one byte per cycle + two output bytes per cycle across segments
    // Approximate as: num_cycles (input) + 2*total_sequence_length (outputs)
    size_t bytes_per_cluster = static_cast<size_t>(num_cycles) + static_cast<size_t>(2 * total_sequence_length);
    // Additional overhead: pointer arrays and read structure
    // Overhead includes pointer arrays and a small fixed margin (~1MB)
    size_t overhead = static_cast<size_t>(num_cycles) * (sizeof(char*) + sizeof(int)) + (1ULL << 20);
    // Use 70% of free memory to be conservative
    size_t usable = static_cast<size_t>(free_mem * 0.70);
    size_t max_clusters_by_mem = usable > overhead && bytes_per_cluster > 0
        ? (usable - overhead) / bytes_per_cluster
        : 0;

    // Determine initial batch size. Prefer memory-derived estimate with sane bounds.
    size_t batch_size = 0;
    if (max_clusters_by_mem > 0) {
        batch_size = std::min<size_t>(num_clusters, max_clusters_by_mem);
    }
    // Lower bound to ensure progress even with tight memory
    const size_t kMinBatch = 16384ULL; // 16k clusters
    if (batch_size == 0) batch_size = std::min<size_t>(num_clusters, kMinBatch);
    // Optional override via env var
    if (const char* env = std::getenv("CUDA_DEMUX_BATCH_SIZE")) {
        try {
            size_t override = std::stoull(env);
            if (override > 0) batch_size = std::min<size_t>(num_clusters, override);
        } catch (...) {}
    }
    std::cout << "Decoding on GPU in batches of up to " << batch_size
              << " clusters (" << num_clusters << " total)." << std::endl;

    // Respect opt-out for adaptive probing via env var
    bool adaptive_probe = true;
    if (const char* no_adapt = std::getenv("CUDA_DEMUX_NO_ADAPTIVE")) {
        if (no_adapt[0] != '\0' && no_adapt[0] != '0') adaptive_probe = false;
    }

    // Adaptive probing: attempt to validate batch size by allocating representative buffers.
    auto can_allocate_for = [&](size_t clusters) -> bool {
        // Inputs: per-cycle this many bytes, outputs: 2 * total_sequence_length bytes per cluster
        size_t in_bytes_total  = static_cast<size_t>(num_cycles) * clusters * sizeof(char);
        size_t out_bytes_total = clusters * static_cast<size_t>(total_sequence_length) * sizeof(char);
        // Temporary buffers for probe
        std::vector<void*> temps;
        temps.reserve(num_cycles + 2);
        // Try allocate outputs first (larger)
        void* tmp_out1 = nullptr; if (cudaMalloc(&tmp_out1, out_bytes_total) != cudaSuccess) { return false; }
        temps.push_back(tmp_out1);
        void* tmp_out2 = nullptr; if (cudaMalloc(&tmp_out2, out_bytes_total) != cudaSuccess) { for (void* p: temps) cudaFree(p); return false; }
        temps.push_back(tmp_out2);
        // Try allocate a single cycle buffer as a proxy for repeated smaller allocations
        void* tmp_in = nullptr; if (cudaMalloc(&tmp_in, clusters * sizeof(char)) != cudaSuccess) { for (void* p: temps) cudaFree(p); return false; }
        temps.push_back(tmp_in);
        // Free probe buffers
        for (void* p: temps) cudaFree(p);
        return true;
    };

    // Binary search down on OOM to find viable batch size, bounded by minimum
    size_t low = std::min<size_t>(batch_size, kMinBatch);
    size_t high = batch_size;
    if (adaptive_probe && !can_allocate_for(batch_size)) {
        std::cerr << "Info: initial batch " << batch_size << " too large; probing smaller sizes." << std::endl;
        while (high > kMinBatch) {
            size_t mid = std::max(kMinBatch, high / 2);
            if (mid == high) break;
            if (can_allocate_for(mid)) { high = mid; break; }
            high = mid;
        }
        // If still not allocable, step down linearly to min
        while (high > kMinBatch && !can_allocate_for(high)) {
            high = std::max(kMinBatch, high / 2);
        }
        batch_size = std::max(kMinBatch, high);
        std::cout << "Selected probed batch size: " << batch_size << " clusters." << std::endl;
    }

    // Ensure inputs make sense
    for (int i = 0; i < num_cycles; ++i) {
        if (h_bcl_sizes[i] < num_clusters) {
            std::cerr << "Warning: cycle " << i << " buffer smaller than num_clusters (" << h_bcl_sizes[i] << " < " << num_clusters << ")" << std::endl;
        }
    }

    // Allocate device-side constant data used across batches
    char** d_bcl_data_ptrs = nullptr; // device array of per-cycle pointers
    int* d_read_structure = nullptr;
    CUDA_CHECK(cudaMalloc(&d_bcl_data_ptrs, num_cycles * sizeof(char*)));
    CUDA_CHECK(cudaMalloc(&d_read_structure, num_cycles * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_read_structure, h_read_structure.data(), num_cycles * sizeof(int), cudaMemcpyHostToDevice));

    // Host output is assembled per batch directly into reads
    reads.resize(num_clusters);

    std::cout << "Decoding on GPU in batches of up to " << batch_size << " clusters (" << num_clusters << " total)." << std::endl;

    // Launch configuration
    int threads_per_block = 256;

    // Buffers allocated per batch
    std::vector<char*> d_bcl_data_buffers(num_cycles, nullptr);

    cudaStream_t stream = nullptr;
    cudaStreamCreate(&stream);

    for (size_t start = 0; start < num_clusters; start += batch_size) {
        size_t this_batch = std::min(batch_size, num_clusters - start);
        
        std::cout << "Processing batch: start=" << start << ", this_batch=" << this_batch << std::endl;

        // In-loop allocation with fallback on OOM by halving this_batch
        char* d_output_sequences = nullptr;
        char* d_output_qualities = nullptr;
        size_t out_bytes = this_batch * static_cast<size_t>(total_sequence_length) * sizeof(char);

        while (true) {
            bool ok = true;
            // Free any previous buffers before retry
            for (int i = 0; i < num_cycles; ++i) { if (d_bcl_data_buffers[i]) { cudaFree(d_bcl_data_buffers[i]); d_bcl_data_buffers[i]=nullptr; } }
            if (d_output_sequences) { cudaFree(d_output_sequences); d_output_sequences = nullptr; }
            if (d_output_qualities) { cudaFree(d_output_qualities); d_output_qualities = nullptr; }

            // Allocate per-cycle input buffers for this batch and copy slices
            for (int i = 0; i < num_cycles; ++i) {
                if (h_bcl_sizes[i] < start + this_batch) {
                    std::cerr << "ERROR: Cycle " << i << " buffer too small: " << h_bcl_sizes[i]
                              << " < " << (start + this_batch) << std::endl;
                    return;
                }
                cudaError_t e = cudaMalloc(&d_bcl_data_buffers[i], this_batch * sizeof(char));
                if (e != cudaSuccess) { ok = false; break; }
            }
            if (ok) {
                // Allocate outputs
                cudaError_t e1 = cudaMalloc(&d_output_sequences, out_bytes);
                cudaError_t e2 = (e1 == cudaSuccess) ? cudaMalloc(&d_output_qualities, out_bytes) : cudaErrorMemoryAllocation;
                ok = (e1 == cudaSuccess && e2 == cudaSuccess);
            }
            if (ok) {
                // Initialize outputs to safe defaults in case some cycles are skipped/missing
                // Sequences -> 'N' (78), Qualities -> '!' (33)
                cudaError_t m1 = cudaMemsetAsync(d_output_sequences, 'N', out_bytes, stream);
                cudaError_t m2 = cudaMemsetAsync(d_output_qualities, '!', out_bytes, stream);
                if (m1 != cudaSuccess || m2 != cudaSuccess) { ok = false; }
            }

            if (ok) {
                // Perform H2D copies (async) and set device pointer array
                for (int i = 0; i < num_cycles; ++i) {
                    const char* h_src = h_bcl_data[i] + start;
                    cudaError_t e = cudaMemcpyAsync(d_bcl_data_buffers[i], h_src, this_batch * sizeof(char), cudaMemcpyHostToDevice, stream);
                    if (e != cudaSuccess) { ok = false; break; }
                }
                if (ok) {
                    cudaError_t e = cudaMemcpyAsync(d_bcl_data_ptrs, d_bcl_data_buffers.data(), num_cycles * sizeof(char*), cudaMemcpyHostToDevice, stream);
                    if (e != cudaSuccess) ok = false;
                }
            }

            if (ok) {
                // Sync to ensure H2D complete before kernel
                cudaError_t e = cudaStreamSynchronize(stream);
                if (e == cudaSuccess) break; else ok = false;
            }

            // On failure, free and halve
            for (int i = 0; i < num_cycles; ++i) { if (d_bcl_data_buffers[i]) { cudaFree(d_bcl_data_buffers[i]); d_bcl_data_buffers[i]=nullptr; } }
            if (d_output_sequences) { cudaFree(d_output_sequences); d_output_sequences = nullptr; }
            if (d_output_qualities) { cudaFree(d_output_qualities); d_output_qualities = nullptr; }

            if (this_batch <= kMinBatch) {
                std::cerr << "ERROR: Unable to allocate for minimum batch size (" << kMinBatch << "). Aborting lane." << std::endl;
                return;
            }
            size_t prev = this_batch;
            this_batch = std::max(kMinBatch, this_batch / 2);
            out_bytes = this_batch * static_cast<size_t>(total_sequence_length) * sizeof(char);
            std::cerr << "Info: Downsizing batch due to allocation failure: " << prev << " -> " << this_batch << std::endl;
        }

        std::cout << "Allocated " << num_cycles << " device buffers, each with " << this_batch << " bytes" << std::endl;

        // Launch kernel
        int blocks_per_grid = static_cast<int>((this_batch + static_cast<size_t>(threads_per_block) - 1) / static_cast<size_t>(threads_per_block));
        std::cout << "Launching kernel: blocks=" << blocks_per_grid 
                  << ", threads_per_block=" << threads_per_block 
                  << ", num_clusters=" << static_cast<int>(this_batch) << std::endl;
        decode_bcl_kernel<<<blocks_per_grid, threads_per_block, 0, stream>>>(
            (const char**)d_bcl_data_ptrs,
            d_read_structure,
            d_output_sequences,
            d_output_qualities,
            num_cycles,
            total_sequence_length,
            static_cast<size_t>(this_batch),
            read1_len,
            index1_len,
            index2_len,
            read2_len
        );
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaStreamSynchronize(stream));

        // Copy results back
        // Use pinned host buffers for faster D2H
        char* h_pinned_seq = nullptr;
        char* h_pinned_qual = nullptr;
        CUDA_CHECK(cudaHostAlloc((void**)&h_pinned_seq, out_bytes, cudaHostAllocDefault));
        CUDA_CHECK(cudaHostAlloc((void**)&h_pinned_qual, out_bytes, cudaHostAllocDefault));
        CUDA_CHECK(cudaMemcpyAsync(h_pinned_seq, d_output_sequences, out_bytes, cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(h_pinned_qual, d_output_qualities, out_bytes, cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));

        // Free batch outputs on device
        cudaFree(d_output_sequences);
        cudaFree(d_output_qualities);

        // Assemble reads for this batch
        for (size_t j = 0; j < this_batch; ++j) {
            size_t read_idx = start + j;
            const char* seq_ptr = h_pinned_seq + j * total_sequence_length;
            const char* qual_ptr = h_pinned_qual + j * total_sequence_length;

            int offset = 0;
            reads[read_idx].sequence.assign(seq_ptr + offset, read1_len);
            reads[read_idx].quality.assign(qual_ptr + offset, read1_len);
            offset += read1_len;

            reads[read_idx].index1.assign(seq_ptr + offset, index1_len);
            offset += index1_len;

            reads[read_idx].index2.assign(seq_ptr + offset, index2_len);
            offset += index2_len;

            reads[read_idx].read2_sequence.assign(seq_ptr + offset, read2_len);
            reads[read_idx].read2_quality.assign(qual_ptr + offset, read2_len);
        }

        // Free per-cycle input buffers for this batch
        for (int i = 0; i < num_cycles; ++i) {
            cudaFree(d_bcl_data_buffers[i]);
            d_bcl_data_buffers[i] = nullptr;
        }
        cudaFreeHost(h_pinned_seq);
        cudaFreeHost(h_pinned_qual);
        std::cout << "  Completed batch " << (start / batch_size + 1) << "/" << ((num_clusters + batch_size - 1) / batch_size)
                  << " (" << this_batch << " clusters)." << std::endl;

        // Gentle self-tuning: re-estimate usable memory and adjust next batch modestly
        size_t free_now=0, total_now=0;
        if (cudaMemGetInfo(&free_now, &total_now) == cudaSuccess && bytes_per_cluster > 0) {
            size_t usable_now = static_cast<size_t>(free_now * 0.70);
            size_t max_by_mem = usable_now > overhead ? (usable_now - overhead) / bytes_per_cluster : batch_size;
            size_t target = std::min(batch_size + batch_size/10, max_by_mem); // +10% up to limit
            if (target >= kMinBatch) batch_size = target;
        }
    }

    // Free persistent device resources
    cudaFree(d_bcl_data_ptrs);
    cudaFree(d_read_structure);
    if (stream) cudaStreamDestroy(stream);

    std::cout << "Assembled " << reads.size() << " reads from GPU results (streamed)." << std::endl;
}
