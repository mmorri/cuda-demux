#include "demux.h"
#include <cuda_runtime.h>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>

struct Read {
    std::string sequence;
    std::string quality;
};

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

std::unordered_map<std::string, std::vector<Read>> demux(const std::vector<Read>& reads, const std::string& samplesheet) {
    std::unordered_map<std::string, std::vector<Read>> demuxed_data;

    auto barcodes = load_barcodes(samplesheet);
    int num_reads = reads.size();
    int num_barcodes = barcodes.size();
    int barcode_length = barcodes[0].length();

    char* d_reads;
    char* d_barcodes;
    int* d_matches;

    cudaMalloc(&d_reads, num_reads * barcode_length * sizeof(char));
    cudaMalloc(&d_barcodes, num_barcodes * barcode_length * sizeof(char));
    cudaMalloc(&d_matches, num_reads * sizeof(int));

    cudaMemcpy(d_reads, reads.data(), num_reads * barcode_length * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_barcodes, barcodes.data(), num_barcodes * barcode_length * sizeof(char), cudaMemcpyHostToDevice);

    int threads_per_block = 256;
    int blocks_per_grid = (num_reads + threads_per_block - 1) / threads_per_block;

    barcode_matching_kernel<<<blocks_per_grid, threads_per_block>>>(d_reads, d_barcodes, d_matches, num_reads, num_barcodes, barcode_length);

    std::vector<int> matches(num_reads);
    cudaMemcpy(matches.data(), d_matches, num_reads * sizeof(int), cudaMemcpyDeviceToHost);

    for (int i = 0; i < num_reads; ++i) {
        if (matches[i] != -1) {
            demuxed_data[barcodes[matches[i]]].push_back(reads[i]);
        }
    }

    cudaFree(d_reads);
    cudaFree(d_barcodes);
    cudaFree(d_matches);

    return demuxed_data;
}