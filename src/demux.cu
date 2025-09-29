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
// Note: Only [BCLConvert_Data] or [Data] sections are parsed for indices.
//       Cloud-specific sections like [Cloud_Data] are ignored for barcodes.
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
                // Only treat BCLConvert_Data or Data as barcode-bearing sections
                in_data_section = (current_section == "BCLConvert_Data" || current_section == "Data");
                // Reset header column indices when entering a new section
                sample_id_col = index_col = index2_col = lane_col = -1;
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
            
            // First non-empty line in data section should be headers
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
    // Env override: allow forcing i5 reverse-complement behavior
    if (const char* env = std::getenv("CUDA_DEMUX_I5_RC")) {
        std::string v(env);
        std::transform(v.begin(), v.end(), v.begin(), ::tolower);
        if (v == "1" || v == "true" || v == "yes") {
            std::cout << "Instrument override via CUDA_DEMUX_I5_RC=1 -> i5 RC = true" << std::endl;
            return true;
        }
        if (v == "0" || v == "false" || v == "no") {
            std::cout << "Instrument override via CUDA_DEMUX_I5_RC=0 -> i5 RC = false" << std::endl;
            return false;
        }
    }
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
    
    // Check CUDA availability and set device if requested via env var
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0) {
        std::cerr << "Error: No CUDA-capable GPU found. " << cudaGetErrorString(err) << std::endl;
        std::cerr << "Please ensure NVIDIA drivers and CUDA are properly installed." << std::endl;
        return demuxed_data;
    }
    // Optional device selection via env var
    if (const char* dev_env = std::getenv("CUDA_DEMUX_DEVICE")) {
        try {
            int dev = std::stoi(dev_env);
            if (dev >= 0 && dev < device_count) {
                cudaSetDevice(dev);
            }
        } catch (...) {}
    }
    
    // Get device properties
    cudaDeviceProp prop;
    int cur_dev = 0; cudaGetDevice(&cur_dev);
    cudaGetDeviceProperties(&prop, cur_dev);
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

    // Optional diagnostic: dump observed barcode counts from decoded reads
    if (const char* dump_env = std::getenv("CUDA_DEMUX_DUMP_BARCODE_COUNTS")) {
        long limit = 200000; // default sample size
        try { if (dump_env[0]) limit = std::stol(dump_env); } catch (...) {}
        if (limit <= 0 || limit > static_cast<long>(reads.size())) limit = static_cast<long>(reads.size());
        std::unordered_map<std::string, long long> cnt_i1, cnt_i2, cnt_i1i2, cnt_i1i2rc;
        for (long i = 0; i < limit; ++i) {
            const auto& r = reads[i];
            if (!r.index1.empty()) cnt_i1[r.index1]++;
            if (!r.index2.empty()) cnt_i2[r.index2]++;
            if (!r.index1.empty() && !r.index2.empty()) {
                cnt_i1i2[r.index1 + r.index2]++;
                cnt_i1i2rc[r.index1 + SampleInfo::reverseComplement(r.index2)]++;
            }
        }
        auto dump_top = [&](const char* title, const std::unordered_map<std::string,long long>& m) {
            std::vector<std::pair<std::string,long long>> v(m.begin(), m.end());
            std::sort(v.begin(), v.end(), [](auto& a, auto& b){ return a.second > b.second; });
            std::cout << "Top " << title << " (first " << limit << " reads):" << std::endl;
            int show = std::min<int>(50, static_cast<int>(v.size()));
            long long total = 0; for (auto& p : v) total += p.second;
            for (int idx = 0; idx < show; ++idx) {
                double pct = total ? (100.0 * v[idx].second / total) : 0.0;
                std::cout << "  " << (idx+1) << ". " << v[idx].first << "\t" << v[idx].second << " (" << pct << "%)" << std::endl;
            }
        };
        dump_top("I1", cnt_i1);
        dump_top("I2", cnt_i2);
        dump_top("I1+I2", cnt_i1i2);
        dump_top("I1+rc(I2)", cnt_i1i2rc);
        // Also print first few SampleSheet barcodes for reference
        std::cout << "SampleSheet combined barcodes (first 20):" << std::endl;
        int printed = 0;
        for (const auto& s : samples) {
            std::cout << "  " << s.sample_id << ": " << s.getCombinedBarcode() << std::endl;
            if (++printed >= 20) break;
        }
    }
    
    // Build initial barcodes and maps; optionally we will try both i5 orientations
    auto build_maps = [&](const std::vector<SampleInfo>& ss,
                          std::unordered_map<std::string, SampleInfo>& map_out,
                          std::vector<std::string>& barcodes_out) {
        map_out.clear();
        for (const auto& s : ss) map_out[s.getCombinedBarcode()] = s;
        barcodes_out = get_combined_barcodes(ss);
    };
    std::unordered_map<std::string, SampleInfo> barcode_to_sample;
    std::vector<std::string> barcodes;
    build_maps(samples, barcode_to_sample, barcodes);
    if (barcodes.empty()) { std::cerr << "Error: No valid barcodes extracted from sample information" << std::endl; return demuxed_data; }

    int num_reads = reads.size();
    int num_barcodes = static_cast<int>(barcodes.size());
    int barcode_length = static_cast<int>(barcodes[0].length());

    std::cout << "Processing " << num_reads << " reads with " << num_barcodes << " barcodes of length " << barcode_length << std::endl;

    // Prepare host data for CUDA processing
    std::vector<char> h_reads(static_cast<size_t>(num_reads) * barcode_length);

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

        // Debug: Print first few reads' indices
        if (i < 5) {
            std::cout << "Debug: Read " << i << " index1='" << reads[i].index1
                      << "' index2='" << reads[i].index2 << "'" << std::endl;
        }

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

    // Allocate device memory for reads once
    char* d_reads;
    err = cudaMalloc(&d_reads, static_cast<size_t>(num_reads) * barcode_length * sizeof(char));
    if (err != cudaSuccess) {
        std::cerr << "CUDA Error (reads): " << cudaGetErrorString(err) << std::endl;
        return demuxed_data;
    }

    // Copy data to device
    err = cudaMemcpy(d_reads, h_reads.data(), static_cast<size_t>(num_reads) * barcode_length * sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "CUDA Memcpy Error (reads): " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_reads);
        return demuxed_data;
    }
    auto run_match = [&](const std::vector<std::string>& bcodes,
                         std::vector<int>& out_matches,
                         std::unordered_map<int,std::string>& out_idx2bc,
                         std::unordered_map<int,std::string>& out_idx2sid) -> bool {
        int nb = static_cast<int>(bcodes.size());
        std::vector<char> h_barcodes(static_cast<size_t>(nb) * barcode_length);
        for (int i = 0; i < nb; ++i) {
            for (int j = 0; j < barcode_length; ++j) {
                h_barcodes[i * barcode_length + j] = bcodes[i][j];
            }
        }
        char* d_barcodes = nullptr; int* d_matches = nullptr;
        cudaError_t e1 = cudaMalloc(&d_barcodes, static_cast<size_t>(nb) * barcode_length * sizeof(char));
        if (e1 != cudaSuccess) { std::cerr << "CUDA Error (barcodes): " << cudaGetErrorString(e1) << std::endl; return false; }
        cudaError_t e2 = cudaMalloc(&d_matches, static_cast<size_t>(num_reads) * sizeof(int));
        if (e2 != cudaSuccess) { cudaFree(d_barcodes); std::cerr << "CUDA Error (matches): " << cudaGetErrorString(e2) << std::endl; return false; }
        cudaError_t e3 = cudaMemcpy(d_barcodes, h_barcodes.data(), static_cast<size_t>(nb) * barcode_length * sizeof(char), cudaMemcpyHostToDevice);
        if (e3 != cudaSuccess) { std::cerr << "CUDA Memcpy Error (barcodes): " << cudaGetErrorString(e3) << std::endl; cudaFree(d_reads); cudaFree(d_barcodes); cudaFree(d_matches); return false; }
        int threads_per_block = 256;
        int blocks_per_grid = (num_reads + threads_per_block - 1) / threads_per_block;
        std::cout << "Launching CUDA kernel with " << blocks_per_grid << " blocks, " << threads_per_block << " threads per block" << std::endl;
        barcode_matching_kernel<<<blocks_per_grid, threads_per_block>>>(d_reads, d_barcodes, d_matches, num_reads, nb, barcode_length);
        cudaError_t ker = cudaGetLastError();
        if (ker != cudaSuccess) { std::cerr << "CUDA Kernel Error: " << cudaGetErrorString(ker) << std::endl; cudaFree(d_barcodes); cudaFree(d_matches); return false; }
        out_matches.assign(num_reads, -1);
        cudaError_t e4 = cudaMemcpy(out_matches.data(), d_matches, static_cast<size_t>(num_reads) * sizeof(int), cudaMemcpyDeviceToHost);
        if (e4 != cudaSuccess) { std::cerr << "CUDA Memcpy Error (matches): " << cudaGetErrorString(e4) << std::endl; cudaFree(d_barcodes); cudaFree(d_matches); return false; }
        // Build index maps
        out_idx2bc.clear(); out_idx2sid.clear();
        for (int i = 0; i < nb; ++i) {
            out_idx2bc[i] = bcodes[i];
            auto it = barcode_to_sample.find(bcodes[i]);
            out_idx2sid[i] = (it != barcode_to_sample.end()) ? it->second.sample_id : std::string();
        }
        cudaFree(d_barcodes); cudaFree(d_matches);
        return true;
    };

    // First attempt: detected orientation
    std::vector<int> matches_a; std::unordered_map<int,std::string> idx2bc_a, idx2sid_a;
    if (!run_match(barcodes, matches_a, idx2bc_a, idx2sid_a)) { cudaFree(d_reads); return demuxed_data; }
    auto score_matches = [&](const std::vector<int>& m,
                             const std::unordered_map<int,std::string>& idx2bc,
                             const std::unordered_map<std::string, SampleInfo>& bc2s) {
        long long matched = 0; std::unordered_map<std::string,int> per_sample;
        for (int i = 0; i < num_reads; ++i) {
            int mi = m[i];
            if (mi >= 0) {
                auto itb = idx2bc.find(mi);
                if (itb != idx2bc.end()) {
                    const auto& bc = itb->second; auto its = bc2s.find(bc);
                    if (its != bc2s.end()) {
                        const SampleInfo& sinfo = its->second;
                        int req_lane = sinfo.lane; int read_lane = reads[i].lane;
                        if (req_lane == 0 || req_lane == read_lane) { matched++; per_sample[sinfo.sample_id]++; }
                    }
                }
            }
        }
        return matched;
    };
    long long score_a = score_matches(matches_a, idx2bc_a, barcode_to_sample);

    // Optional second attempt: opposite i5 orientation if env requests trying both or zero matched
    bool try_both = false; if (const char* tb = std::getenv("CUDA_DEMUX_TRY_BOTH_I5")) { if (tb[0] && tb[0] != '0') try_both = true; }
    std::vector<int> matches_b; std::unordered_map<int,std::string> idx2bc_b, idx2sid_b; std::vector<std::string> barcodes_b; std::unordered_map<std::string, SampleInfo> bc2s_b;
    long long score_b = -1;
    if (try_both || score_a == 0) {
        // Build maps with flipped i5
        std::vector<SampleInfo> samples_flip = samples;
        for (auto& s : samples_flip) s.reverse_complement_i2 = !rc_i5;
        build_maps(samples_flip, bc2s_b, barcodes_b);
        if (!barcodes_b.empty()) {
            // Temporarily swap barcode_to_sample for mapping resolution during this run
            auto saved_map = barcode_to_sample; auto saved_barcodes = barcodes;
            barcode_to_sample = bc2s_b; barcodes = barcodes_b;
            if (!run_match(barcodes_b, matches_b, idx2bc_b, idx2sid_b)) { cudaFree(d_reads); return demuxed_data; }
            score_b = score_matches(matches_b, idx2bc_b, barcode_to_sample);
            // Restore
            barcode_to_sample = saved_map; barcodes = saved_barcodes;
        }
    }

    // Choose better orientation
    const std::vector<int>* matches_ptr = &matches_a;
    const std::unordered_map<int,std::string>* idx2bc_ptr = &idx2bc_a;
    const std::unordered_map<int,std::string>* idx2sid_ptr = &idx2sid_a;
    bool using_flip = false;
    if (score_b > score_a) {
        matches_ptr = &matches_b; idx2bc_ptr = &idx2bc_b; idx2sid_ptr = &idx2sid_b; using_flip = true;
        // Update barcode_to_sample to flipped for grouping
        barcode_to_sample = bc2s_b;
        barcodes = barcodes_b;
        std::cout << "Selected i5 orientation: reverse-complement (score " << score_b << ")" << std::endl;
    } else {
        std::cout << "Selected i5 orientation: " << (rc_i5 ? "reverse-complement" : "forward") << " (score " << score_a << ")" << std::endl;
    }

    // Grouping pass with chosen matches
    std::string undetermined = "undetermined";
    long long matched = 0; long long unmatched = 0;
    std::unordered_map<std::string, int> sample_matches;
    for (int i = 0; i < num_reads; ++i) {
        if (is_adapter_sequence(reads[i].sequence)) {
            demuxed_data["adapter_contamination"].push_back(reads[i]);
            unmatched++;
        } else {
            int mi = (*matches_ptr)[i];
            if (mi >= 0) {
                auto itb = idx2bc_ptr->find(mi);
                if (itb != idx2bc_ptr->end()) {
                    const auto& bc = itb->second; auto its = barcode_to_sample.find(bc);
                    if (its != barcode_to_sample.end()) {
                        const SampleInfo& sinfo = its->second;
                        int req_lane = sinfo.lane; int read_lane = reads[i].lane;
                        if (req_lane == 0 || req_lane == read_lane) {
                            demuxed_data[sinfo.sample_id].push_back(reads[i]); sample_matches[sinfo.sample_id]++; matched++; continue;
                        }
                    }
                }
                demuxed_data[undetermined].push_back(reads[i]); unmatched++;
            } else { demuxed_data[undetermined].push_back(reads[i]); unmatched++; }
        }
    }

    // Report statistics on matched samples
    std::cout << "\nSample matching summary:" << std::endl;
    std::cout << "----------------------" << std::endl;
    for (const auto& sample : samples) {
        int count = sample_matches[sample.sample_id];
        std::cout << sample.sample_id << ": " << count << " reads" << std::endl;
    }
    std::cout << "Demultiplexing complete: " << matched << " reads matched, " << unmatched << " reads unmatched" << std::endl;

    cudaFree(d_reads);
    return demuxed_data;
}
