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
    if (doc.LoadFile(run_info_path.string().c_str()) != XML_SUCCESS) {
        std::cerr << "Error: Could not parse RunInfo.xml." << std::endl;
        return {};
    }

    XMLElement* run_element = doc.FirstChildElement("Run");
    if (!run_element) {
        std::cerr << "Error: <Run> element not found in RunInfo.xml." << std::endl;
        return {};
    }

    XMLElement* reads_element = run_element->FirstChildElement("Reads");
    if (!reads_element) {
        std::cerr << "Error: <Reads> element not found in RunInfo.xml." << std::endl;
        return {};
    }

    RunStructure run_structure;
    int total_cycles = 0;
    for (XMLElement* read_elem = reads_element->FirstChildElement("Read"); read_elem != nullptr; read_elem = read_elem->NextSiblingElement("Read")) {
        int num_cycles = std::stoi(read_elem->Attribute("NumCycles"));
        bool is_indexed = (std::string(read_elem->Attribute("IsIndexedRead")) == "Y");
        int segment_type;
        if (!is_indexed) {
            if (run_structure.read1_cycles == 0) { run_structure.read1_cycles = num_cycles; segment_type = 0; }
            else { run_structure.read2_cycles = num_cycles; segment_type = 3; }
        } else {
            if (run_structure.index1_cycles == 0) { run_structure.index1_cycles = num_cycles; segment_type = 1; }
            else { run_structure.index2_cycles = num_cycles; segment_type = 2; }
        }
        for (int i = 0; i < num_cycles; ++i) { run_structure.read_segments.push_back(segment_type); }
        total_cycles += num_cycles;
    }

    std::cout << "Run Structure: R1:" << run_structure.read1_cycles << " I1:" << run_structure.index1_cycles 
              << " I2:" << run_structure.index2_cycles << " R2:" << run_structure.read2_cycles << std::endl;

    // 2. Detect BCL format and parse accordingly
    fs::path basecalls_dir = bcl_dir / "Data" / "Intensities" / "BaseCalls" / "L001";
    bool has_cbcl = false;
    for(const auto& entry : fs::directory_iterator(basecalls_dir)) {
        if (entry.path().extension() == ".cbcl") {
            has_cbcl = true;
            break;
        }
    }

    if (has_cbcl) {
        std::cout << "CBCL format detected." << std::endl;
        return parse_cbcl(bcl_dir, run_structure, total_cycles);
    } else {
        std::cout << "Legacy BCL (.bcl.gz) format detected." << std::endl;
        return parse_legacy_bcl(bcl_dir, run_structure, total_cycles);
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
    for (int c = 1; c <= total_cycles; ++c) {
        fs::path bcl_file = basecalls_dir / ("C" + std::to_string(c) + ".1") / "L001_1.bcl.gz";
        if (!fs::exists(bcl_file)) { 
             bcl_file = basecalls_dir / ("C" + std::to_string(c) + ".1") / "s_1_1101.bcl.gz";
             if(!fs::exists(bcl_file)) throw std::runtime_error("BCL file not found for cycle " + std::to_string(c));
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

    z_stream strm = {};
    strm.avail_in = block_size;
    strm.next_in = (Bytef*)compressed_data.data();
    
    if (inflateInit2(&strm, 16 + MAX_WBITS) != Z_OK) {
        throw std::runtime_error("Failed to initialize zlib stream");
    }

    std::vector<char> decompressed_data(8192); // Start with a reasonable size
    int ret;
    do {
        strm.avail_out = decompressed_data.size() - strm.total_out;
        strm.next_out = (Bytef*)(decompressed_data.data() + strm.total_out);
        ret = inflate(&strm, Z_NO_FLUSH);
        if (ret != Z_OK && ret != Z_STREAM_END) {
            inflateEnd(&strm);
            throw std::runtime_error("Zlib inflation failed with error code: " + std::to_string(ret));
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

    // 1. Read .filter file to get passing cluster indices
    fs::path filter_file = basecalls_dir / "s_1.filter";
    if (!fs::exists(filter_file)) throw std::runtime_error("Filter file not found: " + filter_file.string());
    
    std::ifstream filter_stream(filter_file, std::ios::binary);
    uint32_t filter_version, num_clusters_total;
    filter_stream.read(reinterpret_cast<char*>(&filter_version), sizeof(filter_version));
    filter_stream.read(reinterpret_cast<char*>(&num_clusters_total), sizeof(num_clusters_total));
    
    std::vector<bool> passing_clusters(num_clusters_total);
    uint32_t num_clusters_passed = 0;
    for (uint32_t i = 0; i < num_clusters_total; ++i) {
        char passed_char;
        filter_stream.read(&passed_char, 1);
        passing_clusters[i] = (passed_char == 1);
        if (passing_clusters[i]) num_clusters_passed++;
    }
    std::cout << "Found " << num_clusters_total << " total clusters, with " << num_clusters_passed << " passing QC." << std::endl;

    // 2. Find and sort all CBCL files
    std::vector<fs::path> cbcl_files;
    for (const auto& entry : fs::directory_iterator(basecalls_dir)) {
        if (entry.path().extension() == ".cbcl") {
            cbcl_files.push_back(entry.path());
        }
    }
    std::sort(cbcl_files.begin(), cbcl_files.end());

    // 3. Prepare buffers for final, cycle-major, filtered BCL data
    std::vector<std::vector<char>> h_bcl_data_buffers(total_cycles, std::vector<char>(num_clusters_passed));
    
    // 4. Process all CBCL files
    uint32_t current_cycle_offset = 0;
    for (const auto& cbcl_path : cbcl_files) {
        std::ifstream cbcl_file(cbcl_path, std::ios::binary);
        if (!cbcl_file) throw std::runtime_error("Failed to open CBCL file: " + cbcl_path.string());

        uint32_t header_size, version, num_gzip_blocks;
        uint32_t bits_per_basecall, bits_per_qscore;
        cbcl_file.read(reinterpret_cast<char*>(&header_size), sizeof(header_size));
        cbcl_file.read(reinterpret_cast<char*>(&version), sizeof(version));
        cbcl_file.seekg(9);
        cbcl_file.read(reinterpret_cast<char*>(&bits_per_basecall), 1);
        cbcl_file.read(reinterpret_cast<char*>(&bits_per_qscore), 1);
        cbcl_file.seekg(header_size - 4);
        cbcl_file.read(reinterpret_cast<char*>(&num_gzip_blocks), sizeof(num_gzip_blocks));

        std::vector<std::pair<uint32_t, uint32_t>> block_offsets_sizes(num_gzip_blocks);
        for (uint32_t i = 0; i < num_gzip_blocks; ++i) {
            cbcl_file.read(reinterpret_cast<char*>(&block_offsets_sizes[i].second), sizeof(uint32_t)); // size
            cbcl_file.seekg(4, std::ios_base::cur); // skip q-score table size
            cbcl_file.read(reinterpret_cast<char*>(&block_offsets_sizes[i].first), sizeof(uint64_t)); // offset
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
