#include "demux.h"
#include <cuda_runtime.h>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

// Function to load barcodes from samplesheet
std::vector<std::string> load_barcodes(const std::string& samplesheet) {
    std::vector<std::string> barcodes;
    std::ifstream file(samplesheet);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open samplesheet file: " << samplesheet << std::endl;
        return barcodes;
    }
    
    std::string line;
    // Skip header line
    std::getline(file, line);
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string sample_id, barcode;
        
        if (std::getline(iss, sample_id, ',') && std::getline(iss, barcode, ',')) {
            barcodes.push_back(barcode);
        }
    }
    
    return barcodes;
}



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
                    // Early termination if we exceed our threshold
                    if (mismatches > 1) {
                        break;
                    }
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

    // Load barcodes from samplesheet
    std::vector<std::string> barcodes = load_barcodes(samplesheet);
    if (barcodes.empty()) {
        std::cerr << "Error: No barcodes found in samplesheet" << std::endl;
        return demuxed_data;
    }

    int num_reads = reads.size();
    int num_barcodes = barcodes.size();
    int barcode_length = barcodes[0].length();

    std::cout << "Processing " << num_reads << " reads with " << num_barcodes << " barcodes of length " << barcode_length << std::endl;

    // Prepare host data for CUDA processing
    std::vector<char> h_reads(num_reads * barcode_length);
    std::vector<char> h_barcodes(num_barcodes * barcode_length);

    // Extract the barcodes from the beginning of each read
    for (int i = 0; i < num_reads; i++) {
        // Ensure we only copy up to barcode_length characters
        int len = std::min(barcode_length, static_cast<int>(reads[i].sequence.length()));
        for (int j = 0; j < len; j++) {
            h_reads[i * barcode_length + j] = reads[i].sequence[j];
        }
        // Pad with 'N' if needed
        for (int j = len; j < barcode_length; j++) {
            h_reads[i * barcode_length + j] = 'N';
        }
    }

    // Copy barcodes to host buffer
    for (int i = 0; i < num_barcodes; i++) {
        for (int j = 0; j < barcode_length; j++) {
            h_barcodes[i * barcode_length + j] = barcodes[i][j];
        }
    }

    // Allocate device memory
    char* d_reads;
    char* d_barcodes;
    int* d_matches;
    cudaError_t err;

    err = cudaMalloc(&d_reads, num_reads * barcode_length * sizeof(char));
    if (err != cudaSuccess) {
        std::cerr << "CUDA Error (reads): " << cudaGetErrorString(err) << std::endl;
        return demuxed_data;
    }

    err = cudaMalloc(&d_barcodes, num_barcodes * barcode_length * sizeof(char));
    if (err != cudaSuccess) {
        cudaFree(d_reads);
        std::cerr << "CUDA Error (barcodes): " << cudaGetErrorString(err) << std::endl;
        return demuxed_data;
    }

    err = cudaMalloc(&d_matches, num_reads * sizeof(int));
    if (err != cudaSuccess) {
        cudaFree(d_reads);
        cudaFree(d_barcodes);
        std::cerr << "CUDA Error (matches): " << cudaGetErrorString(err) << std::endl;
        return demuxed_data;
    }

    // Copy data to device
    err = cudaMemcpy(d_reads, h_reads.data(), num_reads * barcode_length * sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "CUDA Memcpy Error (reads): " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_reads);
        cudaFree(d_barcodes);
        cudaFree(d_matches);
        return demuxed_data;
    }

    err = cudaMemcpy(d_barcodes, h_barcodes.data(), num_barcodes * barcode_length * sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "CUDA Memcpy Error (barcodes): " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_reads);
        cudaFree(d_barcodes);
        cudaFree(d_matches);
        return demuxed_data;
    }

    // Launch kernel
    int threads_per_block = 256;
    int blocks_per_grid = (num_reads + threads_per_block - 1) / threads_per_block;

    std::cout << "Launching CUDA kernel with " << blocks_per_grid << " blocks, " << threads_per_block << " threads per block" << std::endl;
    barcode_matching_kernel<<<blocks_per_grid, threads_per_block>>>(d_reads, d_barcodes, d_matches, num_reads, num_barcodes, barcode_length);

    // Check for kernel errors
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA Kernel Error: " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_reads);
        cudaFree(d_barcodes);
        cudaFree(d_matches);
        return demuxed_data;
    }

    // Copy results back to host
    std::vector<int> matches(num_reads, -1);
    err = cudaMemcpy(matches.data(), d_matches, num_reads * sizeof(int), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        std::cerr << "CUDA Memcpy Error (matches): " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_reads);
        cudaFree(d_barcodes);
        cudaFree(d_matches);
        return demuxed_data;
    }

    // Process the matching results
    std::unordered_map<int, std::string> index_to_barcode;
    for (int i = 0; i < barcodes.size(); ++i) {
        index_to_barcode[i] = barcodes[i];
    }
    
    // Create an "undetermined" category
    std::string undetermined = "undetermined";
    
    // Group the reads by their barcode matches
    int matched = 0;
    int unmatched = 0;
    for (int i = 0; i < num_reads; ++i) {
        if (matches[i] != -1) {
            demuxed_data[index_to_barcode[matches[i]]].push_back(reads[i]);
            matched++;
        } else {
            demuxed_data[undetermined].push_back(reads[i]);
            unmatched++;
        }
    }

    std::cout << "Demultiplexing complete: " << matched << " reads matched, " 
              << unmatched << " reads unmatched" << std::endl;
    
    // Free device memory
    cudaFree(d_reads);
    cudaFree(d_barcodes);
    cudaFree(d_matches);

    return demuxed_data;
}