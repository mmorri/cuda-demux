#include "demux.h"
#include <cuda_runtime.h>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

// Function to load sample information from samplesheet
std::vector<SampleInfo> load_sample_info(const std::string& samplesheet) {
    std::vector<SampleInfo> samples;
    std::ifstream file(samplesheet);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open samplesheet file: " << samplesheet << std::endl;
        return samples;
    }
    
    std::string line;
    // Skip header line
    std::getline(file, line);
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string sample_id, index1, index2;
        
        // Parse the CSV format - expect Sample_ID, Index1, Index2
        if (std::getline(iss, sample_id, ',') && 
            std::getline(iss, index1, ',') && 
            std::getline(iss, index2, ',')) {
            
            // Trim any whitespace
            sample_id.erase(0, sample_id.find_first_not_of(" \t\r\n"));
            sample_id.erase(sample_id.find_last_not_of(" \t\r\n") + 1);
            index1.erase(0, index1.find_first_not_of(" \t\r\n"));
            index1.erase(index1.find_last_not_of(" \t\r\n") + 1);
            index2.erase(0, index2.find_first_not_of(" \t\r\n"));
            index2.erase(index2.find_last_not_of(" \t\r\n") + 1);
            
            // Create the sample info
            SampleInfo sample;
            sample.sample_id = sample_id;
            sample.index1 = index1;
            sample.index2 = index2;
            
            samples.push_back(sample);
        }
    }
    
    std::cout << "Loaded " << samples.size() << " samples from samplesheet" << std::endl;
    for (const auto& sample : samples) {
        std::cout << "Sample: " << sample.sample_id 
                  << ", Index1: " << sample.index1 
                  << ", Index2: " << sample.index2 
                  << ", Combined: " << sample.getCombinedBarcode() << std::endl;
    }
    
    return samples;
}

// Function to extract barcodes from sample info
std::vector<std::string> get_combined_barcodes(const std::vector<SampleInfo>& samples) {
    std::vector<std::string> barcodes;
    
    for (const auto& sample : samples) {
        barcodes.push_back(sample.getCombinedBarcode());
    }
    
    return barcodes;
}

// Function to check if a string is an adapter sequence
bool is_adapter_sequence(const std::string& sequence) {
    for (const auto& adapter : COMMON_ADAPTERS) {
        if (sequence.find(adapter) != std::string::npos) {
            return true;
        }
    }
    return false;
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

    // Load sample information from samplesheet
    std::vector<SampleInfo> samples = load_sample_info(samplesheet);
    if (samples.empty()) {
        std::cerr << "Error: No samples found in samplesheet" << std::endl;
        return demuxed_data;
    }
    
    // Create a lookup map from barcode to sample_id
    std::unordered_map<std::string, std::string> barcode_to_sample;
    for (const auto& sample : samples) {
        barcode_to_sample[sample.getCombinedBarcode()] = sample.sample_id;
    }
    
    // Get the combined barcodes for CUDA processing
    std::vector<std::string> barcodes = get_combined_barcodes(samples);
    if (barcodes.empty()) {
        std::cerr << "Error: No valid barcodes extracted from sample information" << std::endl;
        return demuxed_data;
    }

    int num_reads = reads.size();
    int num_barcodes = barcodes.size();
    int barcode_length = barcodes[0].length();

    std::cout << "Processing " << num_reads << " reads with " << num_barcodes << " barcodes of length " << barcode_length << std::endl;

    // Prepare host data for CUDA processing
    std::vector<char> h_reads(num_reads * barcode_length);
    std::vector<char> h_barcodes(num_barcodes * barcode_length);

    // For dual index processing, we need to extract two barcode regions from each read
    // Typically, Index1 is at the beginning of the read and Index2 may be in a different location
    // or in a separate index read. For this implementation, we'll concatenate the indices.
    
    // First, make sure none of the read sequences match adapter sequences
    int adapter_matches = 0;
    for (int i = 0; i < num_reads; i++) {
        if (is_adapter_sequence(reads[i].sequence)) {
            adapter_matches++;
        }
    }
    
    if (adapter_matches > 0) {
        std::cout << "Warning: " << adapter_matches << " reads contain adapter sequences" << std::endl;
    }
    
    // Now extract the index sequences from each read for matching
    // For this implementation, we'll assume:
    // - The first index (index1) starts at position 0 with a standard length
    // - The second index (index2) follows immediately after index1
    // A real implementation would need to handle more complex scenarios
    
    int index1_length = samples[0].index1.length();
    int index2_length = samples[0].index2.length();
    
    // Make sure the barcode length matches our combined indices
    if (barcode_length != index1_length + index2_length) {
        std::cerr << "Error: Barcode length mismatch. Combined indices are " 
                  << index1_length + index2_length << " bp, but expected " 
                  << barcode_length << " bp" << std::endl;
    }
    
    std::cout << "Processing reads using dual indices: Index1 length=" << index1_length 
              << ", Index2 length=" << index2_length << std::endl;
              
    // Extract the barcodes from each read
    for (int i = 0; i < num_reads; i++) {
        // Extract Index1 from the start of the read
        int pos = 0;
        int len1 = std::min(index1_length, static_cast<int>(reads[i].sequence.length()));
        
        for (int j = 0; j < len1; j++) {
            h_reads[i * barcode_length + j] = reads[i].sequence[j];
            pos++;
        }
        
        // Pad with 'N' if needed
        for (int j = len1; j < index1_length; j++) {
            h_reads[i * barcode_length + j] = 'N';
            pos++;
        }
        
        // Extract Index2 (in a real implementation, this might come from a different location)
        // For now, we'll assume it follows Index1 immediately
        int len2 = std::min(index2_length, static_cast<int>(reads[i].sequence.length() - index1_length));
        
        for (int j = 0; j < len2; j++) {
            h_reads[i * barcode_length + index1_length + j] = reads[i].sequence[index1_length + j];
            pos++;
        }
        
        // Pad with 'N' if needed
        for (int j = len2; j < index2_length; j++) {
            h_reads[i * barcode_length + index1_length + j] = 'N';
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
    std::unordered_map<int, std::string> index_to_sample_id;
    
    for (int i = 0; i < barcodes.size(); ++i) {
        index_to_barcode[i] = barcodes[i];
        // Map the barcode to its sample ID
        index_to_sample_id[i] = barcode_to_sample[barcodes[i]];
    }
    
    // Create an "undetermined" category
    std::string undetermined = "undetermined";
    
    // Group the reads by their sample ID (not by barcode)
    int matched = 0;
    int unmatched = 0;
    
    // Keep track of which samples were matched
    std::unordered_map<std::string, int> sample_matches;
    
    for (int i = 0; i < num_reads; ++i) {
        // Check if this read contains adapter sequences
        if (is_adapter_sequence(reads[i].sequence)) {
            // Put adapter sequences in a special category
            demuxed_data["adapter_contamination"].push_back(reads[i]);
            unmatched++;
        }
        // Valid match to a sample's barcode
        else if (matches[i] != -1) {
            // Use the sample ID instead of the barcode as the key
            std::string sample_id = index_to_sample_id[matches[i]];
            demuxed_data[sample_id].push_back(reads[i]);
            sample_matches[sample_id]++;
            matched++;
        } 
        // No match, put in undetermined
        else {
            demuxed_data[undetermined].push_back(reads[i]);
            unmatched++;
        }
    }
    
    // Report statistics on matched samples
    std::cout << "\nSample matching summary:" << std::endl;
    std::cout << "----------------------" << std::endl;
    for (const auto& sample : samples) {
        int count = sample_matches[sample.sample_id];
        std::cout << sample.sample_id << ": " << count << " reads" << std::endl;
    }

    std::cout << "Demultiplexing complete: " << matched << " reads matched, " 
              << unmatched << " reads unmatched" << std::endl;
    
    // Free device memory
    cudaFree(d_reads);
    cudaFree(d_barcodes);
    cudaFree(d_matches);

    return demuxed_data;
}