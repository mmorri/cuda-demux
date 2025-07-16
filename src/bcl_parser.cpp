#include <algorithm>
#include "bcl_parser.h"
#include "bcl_parser_cuda.h"
#include "common.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <stdexcept>
#include <zlib.h>
#include <tinyxml2.h>
#include <omp.h>

namespace fs = std::filesystem;
using namespace tinyxml2;

// Represents the structure of a sequencing run
struct RunStructure {
    int read1_cycles = 0;
    int index1_cycles = 0;
    int index2_cycles = 0;
    int read2_cycles = 0;
    std::vector<int> read_segments; // 0:R1, 1:I1, 2:I2, 3:R2
};

// Forward declarations for parsing functions
std::vector<Read> parse_legacy_bcl(const fs::path& bcl_dir, const RunStructure& run_structure, int total_cycles);
std::vector<Read> parse_cbcl(const fs::path& bcl_dir, const RunStructure& run_structure, int total_cycles);

// Main C++ function to parse BCL data
std::vector<Read> parse_bcl(const std::string& bcl_folder) {
    fs::path bcl_dir(bcl_folder);
    fs::path run_info_path = bcl_dir / "RunInfo.xml";

    if (!fs::exists(run_info_path)) {
        std::cerr << "Error: RunInfo.xml not found in " << bcl_folder << std::endl;
        return {};
    }

    // 1. Parse RunInfo.xml to get run structure
    XMLDocument doc;
    XMLError result = doc.LoadFile(run_info_path.string().c_str());
    if (result != XML_SUCCESS) {
        std::cerr << "Error: Could not parse RunInfo.xml. Error: " << doc.ErrorIDToName(result) << std::endl;
        std::cerr << "File path: " << run_info_path.string() << std::endl;
        return {};
    }

    XMLElement* run_info_element = doc.FirstChildElement("RunInfo");
    if (!run_info_element) {
        std::cerr << "Error: <RunInfo> element not found in RunInfo.xml." << std::endl;
        std::cerr << "Available root elements:" << std::endl;
        for (XMLElement* child = doc.FirstChildElement(); child != nullptr; child = child->NextSiblingElement()) {
            std::cerr << "  - " << child->Name() << std::endl;
        }
        return {};
    }

    XMLElement* run_element = run_info_element->FirstChildElement("Run");
    if (!run_element) {
        std::cerr << "Error: <Run> element not found in RunInfo.xml." << std::endl;
        std::cerr << "Available child elements of RunInfo:" << std::endl;
        for (XMLElement* child = run_info_element->FirstChildElement(); child != nullptr; child = child->NextSiblingElement()) {
            std::cerr << "  - " << child->Name() << std::endl;
        }
        return {};
    }

    XMLElement* reads_element = run_element->FirstChildElement("Reads");
    if (!reads_element) {
        std::cerr << "Error: <Reads> element not found in RunInfo.xml." << std::endl;
        std::cerr << "Available child elements of Run:" << std::endl;
        for (XMLElement* child = run_element->FirstChildElement(); child != nullptr; child = child->NextSiblingElement()) {
            std::cerr << "  - " << child->Name() << std::endl;
        }
        return {};
    }

    RunStructure run_structure;
    int total_cycles = 0;
    int read_count = 0;
    
    std::cout << "Parsing Read elements..." << std::endl;
    for (XMLElement* read_elem = reads_element->FirstChildElement("Read"); read_elem != nullptr; read_elem = read_elem->NextSiblingElement("Read")) {
        read_count++;
        
        const char* num_cycles_attr = read_elem->Attribute("NumCycles");
        const char* is_indexed_attr = read_elem->Attribute("IsIndexedRead");
        
        if (!num_cycles_attr || !is_indexed_attr) {
            std::cerr << "Error: Missing required attributes in Read element " << read_count << std::endl;
            return {};
        }
        
        int num_cycles = std::stoi(num_cycles_attr);
        bool is_indexed = (std::string(is_indexed_attr) == "Y");
        
        std::cout << "Read " << read_count << ": NumCycles=" << num_cycles << ", IsIndexed=" << (is_indexed ? "Y" : "N") << std::endl;
        
        int segment_type;
        if (!is_indexed) {
            if (run_structure.read1_cycles == 0) { 
                run_structure.read1_cycles = num_cycles; 
                segment_type = 0; 
                std::cout << "  -> Assigned to Read1" << std::endl;
            }
            else { 
                run_structure.read2_cycles = num_cycles; 
                segment_type = 3; 
                std::cout << "  -> Assigned to Read2" << std::endl;
            }
        } else {
            if (run_structure.index1_cycles == 0) { 
                run_structure.index1_cycles = num_cycles; 
                segment_type = 1; 
                std::cout << "  -> Assigned to Index1" << std::endl;
            }
            else { 
                run_structure.index2_cycles = num_cycles; 
                segment_type = 2; 
                std::cout << "  -> Assigned to Index2" << std::endl;
            }
        }
        
        for (int i = 0; i < num_cycles; ++i) { 
            run_structure.read_segments.push_back(segment_type); 
        }
        total_cycles += num_cycles;
    }
    
    std::cout << "Parsed " << read_count << " Read elements" << std::endl;

    std::cout << "Run Structure: R1:" << run_structure.read1_cycles << " I1:" << run_structure.index1_cycles 
              << " I2:" << run_structure.index2_cycles << " R2:" << run_structure.read2_cycles << std::endl;

    // 2. Detect BCL format and parse accordingly
    fs::path basecalls_dir = bcl_dir / "Data" / "Intensities" / "BaseCalls" / "L001";
    
    if (!fs::exists(basecalls_dir)) {
        std::cerr << "Error: BaseCalls directory not found: " << basecalls_dir.string() << std::endl;
        return {};
    }
    
    std::cout << "Scanning directory: " << basecalls_dir.string() << std::endl;
    
    bool has_cbcl = false;
    bool has_legacy_bcl = false;
    std::vector<fs::path> cbcl_files;
    std::vector<fs::path> legacy_bcl_files;
    
    // First, scan the main directory
    for(const auto& entry : fs::directory_iterator(basecalls_dir)) {
        std::cout << "  Found: " << entry.path().filename().string() << std::endl;
        
        if (entry.path().extension() == ".cbcl") {
            has_cbcl = true;
            cbcl_files.push_back(entry.path());
        } else if (entry.path().extension() == ".gz" && entry.path().stem().extension() == ".bcl") {
            has_legacy_bcl = true;
            legacy_bcl_files.push_back(entry.path());
        }
    }
    
    // Then, scan cycle directories recursively for CBCL files
    std::cout << "Scanning cycle directories for CBCL files..." << std::endl;
    for(const auto& entry : fs::directory_iterator(basecalls_dir)) {
        if (entry.is_directory() && entry.path().filename().string().starts_with("C")) {
            std::cout << "  Scanning cycle directory: " << entry.path().filename().string() << std::endl;
            
            for(const auto& subentry : fs::directory_iterator(entry.path())) {
                std::cout << "    Found: " << subentry.path().filename().string() << std::endl;
                
                if (subentry.path().extension() == ".cbcl") {
                    has_cbcl = true;
                    cbcl_files.push_back(subentry.path());
                    std::cout << "      -> Added CBCL file: " << subentry.path().filename().string() << std::endl;
                } else if (subentry.path().extension() == ".gz" && subentry.path().stem().extension() == ".bcl") {
                    has_legacy_bcl = true;
                    legacy_bcl_files.push_back(subentry.path());
                    std::cout << "      -> Added legacy BCL file: " << subentry.path().filename().string() << std::endl;
                }
            }
        }
    }
    
    std::cout << "Found " << cbcl_files.size() << " CBCL files and " << legacy_bcl_files.size() << " legacy BCL files" << std::endl;
    
    if (has_cbcl) {
        std::cout << "CBCL format detected. Using CBCL parser." << std::endl;
        return parse_cbcl(bcl_dir, run_structure, total_cycles);
    } else if (has_legacy_bcl) {
        std::cout << "Legacy BCL (.bcl.gz) format detected. Using legacy BCL parser." << std::endl;
        return parse_legacy_bcl(bcl_dir, run_structure, total_cycles);
    } else {
        std::cerr << "Error: No BCL files found in " << basecalls_dir.string() << std::endl;
        std::cerr << "Expected either .cbcl files or .bcl.gz files" << std::endl;
        return {};
    }
}

// --- Legacy BCL Parser (.bcl.gz) ---
std::vector<char> read_bcl_gz_file(const fs::path& path, uint32_t& cluster_count) {
    gzFile file = gzopen(path.string().c_str(), "rb");
    if (!file) throw std::runtime_error("Could not open gzipped file: " + path.string());
    if (gzread(file, &cluster_count, sizeof(cluster_count)) != sizeof(cluster_count)) {
        gzclose(file);
        throw std::runtime_error("Failed to read cluster count from " + path.string());
    }
    std::vector<char> buffer(cluster_count);
    int bytes_read = gzread(file, buffer.data(), cluster_count);
    gzclose(file);
    if (bytes_read != static_cast<int>(cluster_count)) throw std::runtime_error("Failed to read full BCL data from " + path.string());
    return buffer;
}

std::vector<Read> parse_legacy_bcl(const fs::path& bcl_dir, const RunStructure& run_structure, int total_cycles) {
    std::vector<char*> h_bcl_data_ptrs;
    std::vector<std::vector<char>> h_bcl_data_buffers;
    std::vector<size_t> h_bcl_sizes;
    uint32_t num_clusters = 0;

    fs::path basecalls_dir = bcl_dir / "Data" / "Intensities" / "BaseCalls" / "L001";
    std::cout << "Starting legacy BCL parsing in directory: " << basecalls_dir.string() << std::endl;
    
    for (int c = 1; c <= total_cycles; ++c) {
        fs::path bcl_file = basecalls_dir / ("C" + std::to_string(c) + ".1") / "L001_1.bcl.gz";
        if (!fs::exists(bcl_file)) { 
             bcl_file = basecalls_dir / ("C" + std::to_string(c) + ".1") / "s_1_1101.bcl.gz";
             if(!fs::exists(bcl_file)) {
                 std::cerr << "Error: BCL file not found for cycle " << c << std::endl;
                 std::cerr << "Tried:" << std::endl;
                 std::cerr << "  - " << (basecalls_dir / ("C" + std::to_string(c) + ".1") / "L001_1.bcl.gz").string() << std::endl;
                 std::cerr << "  - " << (basecalls_dir / ("C" + std::to_string(c) + ".1") / "s_1_1101.bcl.gz").string() << std::endl;
                 
                 // List available files in the cycle directory
                 fs::path cycle_dir = basecalls_dir / ("C" + std::to_string(c) + ".1");
                 if (fs::exists(cycle_dir)) {
                     std::cerr << "Available files in cycle directory:" << std::endl;
                     for(const auto& entry : fs::directory_iterator(cycle_dir)) {
                         std::cerr << "  - " << entry.path().filename().string() << std::endl;
                     }
                 } else {
                     std::cerr << "Cycle directory does not exist: " << cycle_dir.string() << std::endl;
                 }
                 throw std::runtime_error("BCL file not found for cycle " + std::to_string(c));
             }
        }
        uint32_t current_clusters = 0;
        h_bcl_data_buffers.push_back(read_bcl_gz_file(bcl_file, current_clusters));
        if (c == 1) num_clusters = current_clusters;
        else if (current_clusters != num_clusters) throw std::runtime_error("Inconsistent cluster count across BCL files.");
    }

    for(auto& buffer : h_bcl_data_buffers) {
        h_bcl_data_ptrs.push_back(buffer.data());
        h_bcl_sizes.push_back(buffer.size());
    }

    std::cout << "Loaded " << total_cycles << " BCL files for " << num_clusters << " clusters." << std::endl;

    std::vector<Read> reads;
    decode_bcl_data_cuda(h_bcl_data_ptrs, h_bcl_sizes, run_structure.read_segments, reads, num_clusters);
    return reads;
}

// --- CBCL Parser (.cbcl) ---

// Helper to decompress a single gzip block from a file stream
std::vector<char> decompress_block(std::ifstream& file, std::streampos block_start, uint32_t block_size) {
    file.seekg(block_start);
    std::vector<char> compressed_data(block_size);
    file.read(compressed_data.data(), block_size);

    // Check if we actually read the expected amount of data
    if (file.gcount() != static_cast<std::streamsize>(block_size)) {
        throw std::runtime_error("Failed to read expected block size. Expected: " + 
                                std::to_string(block_size) + ", Got: " + std::to_string(file.gcount()));
    }

    z_stream strm = {};
    strm.avail_in = block_size;
    strm.next_in = (Bytef*)compressed_data.data();
    
    // Try different zlib window bits for different compression formats
    int window_bits = 16 + MAX_WBITS; // Default: gzip format
    int init_result = inflateInit2(&strm, window_bits);
    
    if (init_result != Z_OK) {
        // Try raw deflate format if gzip fails
        inflateEnd(&strm);
        strm = {};
        strm.avail_in = block_size;
        strm.next_in = (Bytef*)compressed_data.data();
        window_bits = MAX_WBITS; // Raw deflate format
        init_result = inflateInit2(&strm, window_bits);
        
        if (init_result != Z_OK) {
            throw std::runtime_error("Failed to initialize zlib stream. Error: " + std::to_string(init_result));
        }
    }

    std::vector<char> decompressed_data(32768); // Start with a larger size
    int ret;
    do {
        strm.avail_out = decompressed_data.size() - strm.total_out;
        strm.next_out = (Bytef*)(decompressed_data.data() + strm.total_out);
        ret = inflate(&strm, Z_NO_FLUSH);
        
        if (ret != Z_OK && ret != Z_STREAM_END) {
            inflateEnd(&strm);
            std::string error_msg = "Zlib inflation failed with error code: " + std::to_string(ret);
            switch (ret) {
                case Z_BUF_ERROR: error_msg += " (Z_BUF_ERROR - buffer error)"; break;
                case Z_DATA_ERROR: error_msg += " (Z_DATA_ERROR - data error)"; break;
                case Z_MEM_ERROR: error_msg += " (Z_MEM_ERROR - memory error)"; break;
                case Z_STREAM_ERROR: error_msg += " (Z_STREAM_ERROR - stream error)"; break;
                default: error_msg += " (unknown error)"; break;
            }
            throw std::runtime_error(error_msg);
        }
        
        if (strm.avail_out == 0 && ret != Z_STREAM_END) {
            decompressed_data.resize(decompressed_data.size() * 2);
        }
    } while (ret != Z_STREAM_END);

    decompressed_data.resize(strm.total_out);
    inflateEnd(&strm);
    return decompressed_data;
}

std::vector<Read> parse_cbcl(const fs::path& bcl_dir, const RunStructure& run_structure, int total_cycles) {
    fs::path basecalls_dir = bcl_dir / "Data" / "Intensities" / "BaseCalls" / "L001";
    
    std::cout << "Starting CBCL parsing in directory: " << basecalls_dir.string() << std::endl;

    // 1. Read filter files to get passing cluster indices
    std::vector<fs::path> filter_files;
    
    // Look for per-tile filter files (s_1_TTTT.filter format)
    for (const auto& entry : fs::directory_iterator(basecalls_dir)) {
        if (entry.path().extension() == ".filter" && 
            entry.path().filename().string().starts_with("s_1_")) {
            filter_files.push_back(entry.path());
        }
    }
    
    // If no per-tile filter files found, try the old single filter file
    if (filter_files.empty()) {
        fs::path single_filter_file = basecalls_dir / "s_1.filter";
        if (fs::exists(single_filter_file)) {
            filter_files.push_back(single_filter_file);
        }
    }
    
    if (filter_files.empty()) {
        std::cerr << "Error: No filter files found in " << basecalls_dir.string() << std::endl;
        std::cerr << "Available files in directory:" << std::endl;
        for(const auto& entry : fs::directory_iterator(basecalls_dir)) {
            std::cerr << "  - " << entry.path().filename().string() << std::endl;
        }
        throw std::runtime_error("No filter files found in " + basecalls_dir.string());
    }
    
    std::cout << "Found " << filter_files.size() << " filter files:" << std::endl;
    for (const auto& filter_file : filter_files) {
        std::cout << "  - " << filter_file.filename().string() << std::endl;
    }
    
    // Sort filter files to ensure consistent ordering
    std::sort(filter_files.begin(), filter_files.end());
    
    // Read all filter files and combine the data
    std::vector<bool> passing_clusters;
    uint32_t num_clusters_total = 0;
    uint32_t num_clusters_passed = 0;
    
    for (const auto& filter_file : filter_files) {
        std::ifstream filter_stream(filter_file, std::ios::binary);
        if (!filter_stream) {
            throw std::runtime_error("Failed to open filter file: " + filter_file.string());
        }
        
        uint32_t filter_version, tile_clusters;
        filter_stream.read(reinterpret_cast<char*>(&filter_version), sizeof(filter_version));
        filter_stream.read(reinterpret_cast<char*>(&tile_clusters), sizeof(tile_clusters));
        
        std::cout << "  Reading " << filter_file.filename().string() 
                  << ": version=" << filter_version << ", clusters=" << tile_clusters << std::endl;
        
        // Read the filter data for this tile
        for (uint32_t i = 0; i < tile_clusters; ++i) {
            char passed_char;
            filter_stream.read(&passed_char, 1);
            bool passed = (passed_char == 1);
            passing_clusters.push_back(passed);
            if (passed) num_clusters_passed++;
        }
        
        num_clusters_total += tile_clusters;
    }
    
    std::cout << "Combined filter data: " << num_clusters_total << " total clusters, with " 
              << num_clusters_passed << " passing QC." << std::endl;

    // 2. Find and sort all CBCL files (including in cycle directories)
    std::vector<fs::path> cbcl_files;
    
    // First check main directory
    for (const auto& entry : fs::directory_iterator(basecalls_dir)) {
        if (entry.path().extension() == ".cbcl") {
            cbcl_files.push_back(entry.path());
        }
    }
    
    // Then check cycle directories
    for (const auto& entry : fs::directory_iterator(basecalls_dir)) {
        if (entry.is_directory() && entry.path().filename().string().starts_with("C")) {
            for (const auto& subentry : fs::directory_iterator(entry.path())) {
                if (subentry.path().extension() == ".cbcl") {
                    cbcl_files.push_back(subentry.path());
                }
            }
        }
    }
    
    std::sort(cbcl_files.begin(), cbcl_files.end());
    
    std::cout << "Found " << cbcl_files.size() << " CBCL files:" << std::endl;
    for (const auto& cbcl_file : cbcl_files) {
        std::cout << "  - " << cbcl_file.filename().string() << std::endl;
    }
    
    if (cbcl_files.empty()) {
        throw std::runtime_error("No CBCL files found in " + basecalls_dir.string());
    }

    // 3. Prepare buffers for final, cycle-major, filtered BCL data
    std::vector<std::vector<char>> h_bcl_data_buffers(total_cycles, std::vector<char>(num_clusters_passed));
    
    // 4. Process all CBCL files
    uint32_t current_cycle_offset = 0;
    for (const auto& cbcl_path : cbcl_files) {
        std::cout << "Processing CBCL file: " << cbcl_path.filename().string() << std::endl;
        
        std::ifstream cbcl_file(cbcl_path, std::ios::binary);
        if (!cbcl_file) {
            throw std::runtime_error("Failed to open CBCL file: " + cbcl_path.string());
        }

        // Read and validate header
        uint32_t header_size, version, num_gzip_blocks;
        uint32_t bits_per_basecall, bits_per_qscore;
        
        cbcl_file.read(reinterpret_cast<char*>(&header_size), sizeof(header_size));
        if (!cbcl_file) {
            throw std::runtime_error("Failed to read header size from " + cbcl_path.string());
        }
        
        cbcl_file.read(reinterpret_cast<char*>(&version), sizeof(version));
        if (!cbcl_file) {
            throw std::runtime_error("Failed to read version from " + cbcl_path.string());
        }
        
        std::cout << "  Header size: " << header_size << ", Version: " << version << std::endl;
        
        // Skip to bits per basecall/qscore
        cbcl_file.seekg(9, std::ios_base::cur);
        cbcl_file.read(reinterpret_cast<char*>(&bits_per_basecall), 1);
        cbcl_file.read(reinterpret_cast<char*>(&bits_per_qscore), 1);
        
        std::cout << "  Bits per basecall: " << static_cast<int>(bits_per_basecall) 
                  << ", Bits per qscore: " << static_cast<int>(bits_per_qscore) << std::endl;
        
        // Seek to num_gzip_blocks
        cbcl_file.seekg(header_size - 4, std::ios_base::beg);
        cbcl_file.read(reinterpret_cast<char*>(&num_gzip_blocks), sizeof(num_gzip_blocks));
        
        if (!cbcl_file) {
            throw std::runtime_error("Failed to read num_gzip_blocks from " + cbcl_path.string());
        }
        
        std::cout << "  Number of gzip blocks: " << num_gzip_blocks << std::endl;

        std::vector<std::pair<uint32_t, uint32_t>> block_offsets_sizes(num_gzip_blocks);
        for (uint32_t i = 0; i < num_gzip_blocks; ++i) {
            cbcl_file.read(reinterpret_cast<char*>(&block_offsets_sizes[i].second), sizeof(uint32_t)); // size
            if (!cbcl_file) {
                throw std::runtime_error("Failed to read block " + std::to_string(i) + " size from " + cbcl_path.string());
            }
            
            cbcl_file.seekg(4, std::ios_base::cur); // skip q-score table size
            cbcl_file.read(reinterpret_cast<char*>(&block_offsets_sizes[i].first), sizeof(uint64_t)); // offset
            if (!cbcl_file) {
                throw std::runtime_error("Failed to read block " + std::to_string(i) + " offset from " + cbcl_path.string());
            }
            
            std::cout << "    Block " << i << ": size=" << block_offsets_sizes[i].second 
                      << ", offset=" << block_offsets_sizes[i].first << std::endl;
        }

        std::vector<std::vector<char>> decompressed_blocks(num_gzip_blocks);
        #pragma omp parallel for
        for (uint32_t i = 0; i < num_gzip_blocks; ++i) {
            decompressed_blocks[i] = decompress_block(cbcl_file, block_offsets_sizes[i].first, block_offsets_sizes[i].second);
        }

        // Transpose from cluster-major to cycle-major and filter
        uint32_t cycles_in_this_cbcl = decompressed_blocks[0].size();
        for (uint32_t c = 0; c < cycles_in_this_cbcl; ++c) {
            uint32_t passed_cluster_idx = 0;
            for (uint32_t i = 0; i < num_clusters_total; ++i) {
                if (passing_clusters[i]) {
                    uint32_t block_idx = i / (1000 * 256);
                    uint32_t tile_idx = (i / 256) % 1000;
                    uint32_t cluster_in_tile_idx = i % 256;
                    uint32_t offset = (tile_idx * 256 + cluster_in_tile_idx) * cycles_in_this_cbcl + c;
                    h_bcl_data_buffers[current_cycle_offset + c][passed_cluster_idx++] = decompressed_blocks[block_idx][offset];
                }
            }
        }
        current_cycle_offset += cycles_in_this_cbcl;
    }

    // 5. Prepare data for CUDA kernel
    std::vector<char*> h_bcl_data_ptrs;
    std::vector<size_t> h_bcl_sizes;
    for(auto& buffer : h_bcl_data_buffers) {
        h_bcl_data_ptrs.push_back(buffer.data());
        h_bcl_sizes.push_back(buffer.size());
    }

    std::cout << "Loaded and processed " << total_cycles << " cycles for " << num_clusters_passed << " passed clusters." << std::endl;

    // 6. Call the CUDA function to decode the BCL data
    std::vector<Read> reads;
    decode_bcl_data_cuda(h_bcl_data_ptrs, h_bcl_sizes, run_structure.read_segments, reads, num_clusters_passed);
    return reads;
}
