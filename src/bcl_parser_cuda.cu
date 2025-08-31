#include "bcl_parser_cuda.h"
#include <cuda_runtime.h>
#include <iostream>
#include <numeric>

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
    int num_clusters              // Total number of clusters (reads)
) {
    int cluster_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cluster_idx >= num_clusters) {
        return;
    }

    static const char bases[] = {'A', 'C', 'G', 'T', 'N'};

    for (int cycle = 0; cycle < num_cycles; ++cycle) {
        // Each BCL file buffer contains data for all clusters for one cycle
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
            int output_idx = cluster_idx * total_sequence_length + cycle;
            d_output_sequences[output_idx] = base;
            d_output_qualities[output_idx] = quality_char;
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
    int total_sequence_length = h_read_structure.size();

    // --- Streaming approach: process clusters in GPU-sized batches ---
    // Calculate a safe batch size based on available GPU memory
    size_t free_mem = 0, total_mem = 0;
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));

    // Bytes per cluster on device: one byte per cycle + two output bytes per cycle across segments
    // Approximate as: num_cycles (input) + 2*total_sequence_length (outputs)
    size_t bytes_per_cluster = static_cast<size_t>(num_cycles) + static_cast<size_t>(2 * total_sequence_length);
    // Additional overhead: pointer arrays and read structure
    size_t overhead = static_cast<size_t>(num_cycles) * (sizeof(char*) + sizeof(int)) + 1 << 20; // ~1MB margin
    // Use 70% of free memory to be conservative
    size_t usable = static_cast<size_t>(free_mem * 0.70);
    size_t max_clusters_by_mem = usable > overhead && bytes_per_cluster > 0
        ? (usable - overhead) / bytes_per_cluster
        : 0;

    // Fallback minimum batch size if estimation is too small
    size_t batch_size = std::min<size_t>(num_clusters, std::max<size_t>(1, max_clusters_by_mem));

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

    int read1_len = 0, index1_len = 0, index2_len = 0, read2_len = 0;
    for (int seg : h_read_structure) {
        if (seg == 0) read1_len++;
        else if (seg == 1) index1_len++;
        else if (seg == 2) index2_len++;
        else if (seg == 3) read2_len++;
    }

    std::cout << "Decoding on GPU in batches of up to " << batch_size << " clusters (" << num_clusters << " total)." << std::endl;

    // Launch configuration
    int threads_per_block = 256;

    // Buffers allocated per batch
    std::vector<char*> d_bcl_data_buffers(num_cycles, nullptr);

    for (size_t start = 0; start < num_clusters; start += batch_size) {
        size_t this_batch = std::min(batch_size, num_clusters - start);

        // Allocate per-cycle input buffers for this batch and copy slices
        for (int i = 0; i < num_cycles; ++i) {
            CUDA_CHECK(cudaMalloc(&d_bcl_data_buffers[i], this_batch * sizeof(char)));
            const char* h_src = h_bcl_data[i] + start;
            CUDA_CHECK(cudaMemcpy(d_bcl_data_buffers[i], h_src, this_batch * sizeof(char), cudaMemcpyHostToDevice));
        }

        // Update device pointer array
        CUDA_CHECK(cudaMemcpy(d_bcl_data_ptrs, d_bcl_data_buffers.data(), num_cycles * sizeof(char*), cudaMemcpyHostToDevice));

        // Allocate outputs for this batch
        char* d_output_sequences = nullptr;
        char* d_output_qualities = nullptr;
        size_t out_bytes = this_batch * static_cast<size_t>(total_sequence_length) * sizeof(char);
        CUDA_CHECK(cudaMalloc(&d_output_sequences, out_bytes));
        CUDA_CHECK(cudaMalloc(&d_output_qualities, out_bytes));

        // Launch kernel
        int blocks_per_grid = static_cast<int>((this_batch + threads_per_block - 1) / threads_per_block);
        decode_bcl_kernel<<<blocks_per_grid, threads_per_block>>>(
            (const char**)d_bcl_data_ptrs,
            d_read_structure,
            d_output_sequences,
            d_output_qualities,
            num_cycles,
            total_sequence_length,
            static_cast<int>(this_batch)
        );
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        // Copy results back
        std::vector<char> h_output_sequences_batch(out_bytes);
        std::vector<char> h_output_qualities_batch(out_bytes);
        CUDA_CHECK(cudaMemcpy(h_output_sequences_batch.data(), d_output_sequences, out_bytes, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_output_qualities_batch.data(), d_output_qualities, out_bytes, cudaMemcpyDeviceToHost));

        // Free batch outputs on device
        cudaFree(d_output_sequences);
        cudaFree(d_output_qualities);

        // Assemble reads for this batch
        for (size_t j = 0; j < this_batch; ++j) {
            size_t read_idx = start + j;
            const char* seq_ptr = h_output_sequences_batch.data() + j * total_sequence_length;
            const char* qual_ptr = h_output_qualities_batch.data() + j * total_sequence_length;

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

        std::cout << "  Completed batch " << (start / batch_size + 1) << "/" << ((num_clusters + batch_size - 1) / batch_size)
                  << " (" << this_batch << " clusters)." << std::endl;
    }

    // Free persistent device resources
    cudaFree(d_bcl_data_ptrs);
    cudaFree(d_read_structure);

    std::cout << "Assembled " << reads.size() << " reads from GPU results (streamed)." << std::endl;
}
