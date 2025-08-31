#include "demux.h"
#include <cuda_runtime.h>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <tinyxml2.h>
#include <filesystem>

// Helper function to trim whitespace
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

// Function to load sample information from Illumina samplesheet
// Supports optional Lane column. If lane is empty/0, applies to all lanes.
std::vector<SampleInfo> load_sample_info(const std::string& samplesheet) {
    std::vector<SampleInfo> samples;
    std::ifstream file(samplesheet);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open samplesheet file: " << samplesheet << std::endl;
        return samples;
    }
    
    std::string line;
    std::string current_section;
    bool in_data_section = false;
    int sample_id_col = -1, index_col = -1, index2_col = -1, lane_col = -1;
    
    while (std::getline(file, line)) {
        // Remove any BOM or special characters
        if (!line.empty() && line[0] == '\xEF') {
            line = line.substr(3);
        }
        
        line = trim(line);
        if (line.empty()) continue;
        
        // Check for section headers
        if (line[0] == '[') {
            size_t end = line.find(']');
            if (end != std::string::npos) {
                current_section = line.substr(1, end - 1);
                in_data_section = (current_section == "BCLConvert_Data" || 
                                 current_section == "Data" || 
                                 current_section == "Cloud_Data");
                continue;
            }
        }
        
        // Process data section
        if (in_data_section) {
            std::vector<std::string> fields;
            std::stringstream ss(line);
            std::string field;
            
            while (std::getline(ss, field, ',')) {
                fields.push_back(trim(field));
            }
            
            // First line in data section should be headers
            if (sample_id_col == -1) {
                for (size_t i = 0; i < fields.size(); ++i) {
                    std::string header = fields[i];
                    // Convert to lowercase for comparison
                    std::transform(header.begin(), header.end(), header.begin(), ::tolower);
                    
                    if (header == "sample_id" || header == "sampleid" || header == "sample") {
                        sample_id_col = i;
                    } else if (header == "index" || header == "index1" || header == "i7_index_id") {
                        index_col = i;
                    } else if (header == "index2" || header == "i5_index_id") {
                        index2_col = i;
                    } else if (header == "lane") {
                        lane_col = i;
                    }
                }
                continue;
            }
            
            // Parse sample data
            if (sample_id_col >= 0 && sample_id_col < fields.size()) {
                SampleInfo sample;
                sample.sample_id = fields[sample_id_col];
                
                // Get Index1 if present
                if (index_col >= 0 && index_col < fields.size()) {
                    sample.index1 = fields[index_col];
                }
                
                // Get Index2 if present (might not exist for single-index)
                if (index2_col >= 0 && index2_col < fields.size()) {
                    sample.index2 = fields[index2_col];
                }

                // Lane if present
                if (lane_col >= 0 && lane_col < fields.size()) {
                    std::string lane_str = fields[lane_col];
                    if (!lane_str.empty()) {
                        try {
                            sample.lane = std::stoi(lane_str);
                        } catch (...) { sample.lane = 0; }
                    }
                }
                
                // Skip empty sample IDs or header-like rows
                if (!sample.sample_id.empty() && 
                    sample.sample_id != "Sample_ID" && 
                    sample.sample_id != "SampleID") {
                    samples.push_back(sample);
                }
            }
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

// Simple parser for RunParameters.xml to determine if i5 should be reverse-complemented
static bool detect_reverse_complement_i5(const std::string& run_folder) {
    namespace fs = std::filesystem;
    fs::path rp = fs::path(run_folder) / "RunParameters.xml";
    if (!fs::exists(rp)) {
        std::cerr << "Warning: RunParameters.xml not found at " << rp.string() << ". Assuming no i5 reverse-complement." << std::endl;
        return false;
    }
    tinyxml2::XMLDocument doc;
    if (doc.LoadFile(rp.string().c_str()) != tinyxml2::XML_SUCCESS) {
        std::cerr << "Warning: Failed to parse RunParameters.xml. Assuming no i5 reverse-complement." << std::endl;
        return false;
    }
    std::string text;
    // Try a few known nodes
    const char* keys[] = {"InstrumentName", "InstrumentType", "ApplicationName", "Platform"};
    for (const char* k : keys) {
        auto* el = doc.FirstChildElement(k);
        if (el && el->GetText()) { text = el->GetText(); break; }
    }
    // Fallback: scan document text
    if (text.empty()) {
        tinyxml2::XMLPrinter pr;
        doc.Print(&pr);
        text = pr.CStr();
    }
    // Heuristic: NextSeq, MiniSeq, NovaSeq require i5 reverse complement
    auto lower = text;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    bool rc = (lower.find("nextseq") != std::string::npos) ||
              (lower.find("miniseq") != std::string::npos) ||
              (lower.find("novaseq") != std::string::npos);
    std::cout << "Instrument detection from RunParameters: '" << text << "' -> i5 RC = " << (rc ? "true" : "false") << std::endl;
    return rc;
}

std::unordered_map<std::string, std::vector<Read>> demux(const std::vector<Read>& reads, const std::string& samplesheet, const std::string& run_folder) {
    std::unordered_map<std::string, std::vector<Read>> demuxed_data;
    
    // Check CUDA availability
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0) {
        std::cerr << "Error: No CUDA-capable GPU found. " << cudaGetErrorString(err) << std::endl;
        std::cerr << "Please ensure NVIDIA drivers and CUDA are properly installed." << std::endl;
        return demuxed_data;
    }
    
    // Get device properties
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Using GPU: " << prop.name << " with compute capability " 
              << prop.major << "." << prop.minor << std::endl;

    // Load sample information from samplesheet
    std::vector<SampleInfo> samples = load_sample_info(samplesheet);
    // Set platform i5 reverse-complement on all samples based on RunParameters
    bool rc_i5 = detect_reverse_complement_i5(run_folder);
    for (auto& s : samples) s.reverse_complement_i2 = rc_i5;
    if (samples.empty()) {
        std::cerr << "Error: No samples found in samplesheet" << std::endl;
        return demuxed_data;
    }
    
    // Create a lookup map from barcode to sample_id with lane restriction
    // We keep barcodes vector for CUDA matching, then validate lane post-match
    std::unordered_map<std::string, SampleInfo> barcode_to_sample;
    for (const auto& sample : samples) {
        barcode_to_sample[sample.getCombinedBarcode()] = sample;
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
              
    // Handle single-index case
    bool is_dual_index = !samples[0].index2.empty();
    
    // Extract the barcodes from each read
    for (int i = 0; i < num_reads; i++) {
        // Extract Index1 from the read's index1 field
        int len1 = std::min(index1_length, static_cast<int>(reads[i].index1.length()));
        
        for (int j = 0; j < len1; j++) {
            h_reads[i * barcode_length + j] = reads[i].index1[j];
        }
        
        // Pad with 'N' if needed
        for (int j = len1; j < index1_length; j++) {
            h_reads[i * barcode_length + j] = 'N';
        }
        
        if (is_dual_index) {
            // Extract Index2 from the read's index2 field
            int len2 = std::min(index2_length, static_cast<int>(reads[i].index2.length()));
            
            for (int j = 0; j < len2; j++) {
                h_reads[i * barcode_length + index1_length + j] = reads[i].index2[j];
            }
            
            // Pad with 'N' if needed
            for (int j = len2; j < index2_length; j++) {
                h_reads[i * barcode_length + index1_length + j] = 'N';
            }
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
        auto it = barcode_to_sample.find(barcodes[i]);
        if (it != barcode_to_sample.end()) {
            index_to_sample_id[i] = it->second.sample_id;
        } else {
            index_to_sample_id[i] = std::string();
        }
    }
    
    // Create an "undetermined" category
    std::string undetermined = "undetermined";
    
    // Group the reads by their sample ID, respecting lane restrictions
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
            // Use the sample ID, validating lane restriction if present
            std::string sample_id = index_to_sample_id[matches[i]];
            const SampleInfo& sinfo = barcode_to_sample[index_to_barcode[matches[i]]];
            int required_lane = sinfo.lane; // 0 means all
            int read_lane = reads[i].lane;
            if (required_lane == 0 || required_lane == read_lane) {
                demuxed_data[sample_id].push_back(reads[i]);
                sample_matches[sample_id]++;
                matched++;
            } else {
                demuxed_data[undetermined].push_back(reads[i]);
                unmatched++;
            }
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
