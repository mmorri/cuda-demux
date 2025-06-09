#include "demux.h"
#include <cuda_runtime.h>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>
#include "read_types.h"

/**
 * @brief CUDA kernel for barcode matching.
 */
__global__ void barcode_matching_kernel(const char* reads, const char* barcodes, int* matches, int num_reads, int num_barcodes, int barcode_length) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < num_reads) {
        int best_match = -1;
        int min_mismatches = barcode_length + 1;

        for (int i = 0; i < num_barcodes; ++i) {
            int mismatches = 0;
            for (int j = 0; j < barcode_length; ++j) {
                if (reads[idx * barcode_length + j] != barcodes[i * barcode_length + j]) {
                    mismatches++;
                }
            }
            if (mismatches < min_mismatches) {
                min_mismatches = mismatches;
                best_match = i;
            }
        }

        matches[idx] = (min_mismatches <= 1) ? best_match : -1;
    }
}

/**
 * @brief Demultiplex reads using CUDA-accelerated barcode matching.
 * @param reads Vector of sequencing reads.
 * @param samplesheet Path to the sample sheet CSV.
 * @return Map from barcode/sample to vector of reads.
 * @throws std::runtime_error on CUDA or file errors.
 */
std::unordered_map<std::string, std::vector<Read>> demux(const std::vector<Read>& reads, const std::string& samplesheet) {
    std::unordered_map<std::string, std::vector<Read>> demuxed_data;

    auto barcodes = load_barcodes(samplesheet);
    int num_reads = reads.size();
    int num_barcodes = barcodes.size();
    int barcode_length = barcodes[0].length();

    // Helper for CUDA error checking
    auto cudaCheck = [](cudaError_t result, const char* msg) {
        if (result != cudaSuccess) {
            throw std::runtime_error(std::string("CUDA Error: ") + msg + ": " + cudaGetErrorString(result));
        }
    };

    char* d_reads = nullptr;
    char* d_barcodes = nullptr;
    int* d_matches = nullptr;

    try {
        cudaCheck(cudaMalloc(&d_reads, num_reads * barcode_length * sizeof(char)), "cudaMalloc d_reads");
        cudaCheck(cudaMalloc(&d_barcodes, num_barcodes * barcode_length * sizeof(char)), "cudaMalloc d_barcodes");
        cudaCheck(cudaMalloc(&d_matches, num_reads * sizeof(int)), "cudaMalloc d_matches");

        // Prepare host data in contiguous memory
        std::vector<char> host_reads(num_reads * barcode_length);
        for (int i = 0; i < num_reads; ++i) {
            memcpy(&host_reads[i * barcode_length], reads[i].sequence.data(), barcode_length);
        }
        std::vector<char> host_barcodes(num_barcodes * barcode_length);
        for (int i = 0; i < num_barcodes; ++i) {
            memcpy(&host_barcodes[i * barcode_length], barcodes[i].data(), barcode_length);
        }

        cudaCheck(cudaMemcpy(d_reads, host_reads.data(), num_reads * barcode_length * sizeof(char), cudaMemcpyHostToDevice), "cudaMemcpy d_reads");
        cudaCheck(cudaMemcpy(d_barcodes, host_barcodes.data(), num_barcodes * barcode_length * sizeof(char), cudaMemcpyHostToDevice), "cudaMemcpy d_barcodes");

        int threads_per_block = 256;
        int blocks_per_grid = (num_reads + threads_per_block - 1) / threads_per_block;

        barcode_matching_kernel<<<blocks_per_grid, threads_per_block>>>(d_reads, d_barcodes, d_matches, num_reads, num_barcodes, barcode_length);
        cudaCheck(cudaGetLastError(), "Kernel launch");
        cudaCheck(cudaDeviceSynchronize(), "Kernel sync");

        std::vector<int> matches(num_reads);
        cudaCheck(cudaMemcpy(matches.data(), d_matches, num_reads * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy d_matches");

        for (int i = 0; i < num_reads; ++i) {
            if (matches[i] != -1) {
                demuxed_data[barcodes[matches[i]]].push_back(reads[i]);
            }
        }
    } catch (...) {
        if (d_reads) cudaFree(d_reads);
        if (d_barcodes) cudaFree(d_barcodes);
        if (d_matches) cudaFree(d_matches);
        throw;
    }
    if (d_reads) cudaFree(d_reads);
    if (d_barcodes) cudaFree(d_barcodes);
    if (d_matches) cudaFree(d_matches);

    return demuxed_data;
}