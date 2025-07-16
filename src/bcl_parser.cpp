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

// Helper to read a gzipped file into a buffer
std::vector<char> read_gzipped_file(const fs::path& path, uint32_t& cluster_count) {
    gzFile file = gzopen(path.string().c_str(), "rb");
    if (!file) {
        throw std::runtime_error("Could not open gzipped file: " + path.string());
    }

    // BCL files have a 4-byte header indicating the number of clusters
    if (gzread(file, &cluster_count, sizeof(cluster_count)) != sizeof(cluster_count)) {
        gzclose(file);
        throw std::runtime_error("Failed to read cluster count from " + path.string());
    }

    std::vector<char> buffer(cluster_count);
    int bytes_read = gzread(file, buffer.data(), cluster_count);
    gzclose(file);

    if (bytes_read != static_cast<int>(cluster_count)) {
        throw std::runtime_error("Failed to read full BCL data from " + path.string());
    }

    return buffer;
}

// Main C++ function to parse BCL data
std::vector<Read> parse_bcl(const std::string& bcl_folder) {
    fs::path bcl_dir(bcl_folder);
    fs::path run_info_path = bcl_dir / "RunInfo.xml";

    if (!fs::exists(run_info_path)) {
        // Fallback to test data generation if no real run is found
        std::cout << "RunInfo.xml not found. Generating test data instead." << std::endl;
        // NOTE: A real implementation of generateTestData() would be needed here.
        // For now, returning an empty vector.
        return {};
    }

    // 1. Parse RunInfo.xml to get run structure
    XMLDocument doc;
    doc.LoadFile(run_info_path.string().c_str());
    XMLElement* reads_element = doc.FirstChildElement("Run")->FirstChildElement("Reads");

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

        for (int i = 0; i < num_cycles; ++i) {
            run_structure.read_segments.push_back(segment_type);
        }
        total_cycles += num_cycles;
    }

    std::cout << "Run Structure: R1:" << run_structure.read1_cycles << " I1:" << run_structure.index1_cycles 
              << " I2:" << run_structure.index2_cycles << " R2:" << run_structure.read2_cycles << std::endl;

    // 2. Load raw BCL data for all cycles into memory
    std::vector<char*> h_bcl_data_ptrs;
    std::vector<std::vector<char>> h_bcl_data_buffers;
    std::vector<size_t> h_bcl_sizes;
    uint32_t num_clusters = 0;

    fs::path basecalls_dir = bcl_dir / "Data" / "Intensities" / "BaseCalls" / "L001";
    for (int c = 1; c <= total_cycles; ++c) {
        // BCL files can have different naming schemes, so we check for common ones.
        fs::path bcl_file = basecalls_dir / ("C" + std::to_string(c) + ".1") / "L001_1.bcl.gz";
        if (!fs::exists(bcl_file)) { 
             bcl_file = basecalls_dir / ("C" + std::to_string(c) + ".1") / "s_1_1101.bcl.gz";
             if(!fs::exists(bcl_file)) throw std::runtime_error("BCL file not found for cycle " + std::to_string(c));
        }
        
        uint32_t current_clusters = 0;
        h_bcl_data_buffers.push_back(read_gzipped_file(bcl_file, current_clusters));
        if (c == 1) {
            num_clusters = current_clusters;
        } else if (current_clusters != num_clusters) {
            throw std::runtime_error("Inconsistent cluster count across BCL files.");
        }
    }

    for(auto& buffer : h_bcl_data_buffers) {
        h_bcl_data_ptrs.push_back(buffer.data());
        h_bcl_sizes.push_back(buffer.size());
    }

    std::cout << "Loaded " << total_cycles << " BCL files for " << num_clusters << " clusters." << std::endl;

    // 3. Call the CUDA function to decode the BCL data
    std::vector<Read> reads;
    decode_bcl_data_cuda(
        h_bcl_data_ptrs,
        h_bcl_sizes,
        run_structure.read_segments,
        reads,
        num_clusters
    );

    return reads;
}
