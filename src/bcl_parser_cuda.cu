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

    // --- 1. Allocate GPU Memory ---
    char** d_bcl_data_ptrs; // An array of pointers on the device
    char** d_bcl_data_buffers; // The actual BCL data buffers on the device
    int* d_read_structure;
    char* d_output_sequences;
    char* d_output_qualities;

    CUDA_CHECK(cudaMalloc(&d_bcl_data_ptrs, num_cycles * sizeof(char*)));
    d_bcl_data_buffers = new char*[num_cycles];

    for (int i = 0; i < num_cycles; ++i) {
        CUDA_CHECK(cudaMalloc(&d_bcl_data_buffers[i], h_bcl_sizes[i]));
    }
    CUDA_CHECK(cudaMalloc(&d_read_structure, num_cycles * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_output_sequences, num_clusters * total_sequence_length * sizeof(char)));
    CUDA_CHECK(cudaMalloc(&d_output_qualities, num_clusters * total_sequence_length * sizeof(char)));

    // --- 2. Transfer Data from Host to GPU ---
    for (int i = 0; i < num_cycles; ++i) {
        CUDA_CHECK(cudaMemcpy(d_bcl_data_buffers[i], h_bcl_data[i], h_bcl_sizes[i], cudaMemcpyHostToDevice));
    }
    CUDA_CHECK(cudaMemcpy(d_bcl_data_ptrs, d_bcl_data_buffers, num_cycles * sizeof(char*), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_read_structure, h_read_structure.data(), num_cycles * sizeof(int), cudaMemcpyHostToDevice));

    // --- 3. Launch Kernel ---
    int threads_per_block = 256;
    int blocks_per_grid = (num_clusters + threads_per_block - 1) / threads_per_block;
    
    std::cout << "Launching BCL decode kernel with " << blocks_per_grid << " blocks and " << threads_per_block << " threads." << std::endl;
    decode_bcl_kernel<<<blocks_per_grid, threads_per_block>>>(
        (const char**)d_bcl_data_ptrs,
        d_read_structure,
        d_output_sequences,
        d_output_qualities,
        num_cycles,
        total_sequence_length,
        num_clusters
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // --- 4. Transfer Results from GPU to Host ---
    std::vector<char> h_output_sequences(num_clusters * total_sequence_length);
    std::vector<char> h_output_qualities(num_clusters * total_sequence_length);
    CUDA_CHECK(cudaMemcpy(h_output_sequences.data(), d_output_sequences, h_output_sequences.size() * sizeof(char), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_output_qualities.data(), d_output_qualities, h_output_qualities.size() * sizeof(char), cudaMemcpyDeviceToHost));

    // --- 5. Free GPU Memory ---
    for (int i = 0; i < num_cycles; ++i) {
        cudaFree(d_bcl_data_buffers[i]);
    }
    delete[] d_bcl_data_buffers;
    cudaFree(d_bcl_data_ptrs);
    cudaFree(d_read_structure);
    cudaFree(d_output_sequences);
    cudaFree(d_output_qualities);

    // --- 6. Assemble Reads on Host ---
    reads.resize(num_clusters);
    int read1_len = 0, index1_len = 0, index2_len = 0, read2_len = 0;
    for(int seg : h_read_structure) {
        if(seg == 0) read1_len++;
        else if(seg == 1) index1_len++;
        else if(seg == 2) index2_len++;
        else if(seg == 3) read2_len++;
    }

    for (size_t i = 0; i < num_clusters; ++i) {
        const char* seq_ptr = h_output_sequences.data() + i * total_sequence_length;
        const char* qual_ptr = h_output_qualities.data() + i * total_sequence_length;

        int offset = 0;
        reads[i].sequence.assign(seq_ptr + offset, read1_len);
        reads[i].quality.assign(qual_ptr + offset, read1_len);
        offset += read1_len;

        reads[i].index1.assign(seq_ptr + offset, index1_len);
        offset += index1_len;

        reads[i].index2.assign(seq_ptr + offset, index2_len);
        offset += index2_len;

        reads[i].read2_sequence.assign(seq_ptr + offset, read2_len);
        reads[i].read2_quality.assign(qual_ptr + offset, read2_len);
    }
    std::cout << "Assembled " << reads.size() << " reads from GPU results." << std::endl;
}
