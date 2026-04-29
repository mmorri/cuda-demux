#include <algorithm>
#include "bcl_parser.h"
#include "common.h"
#include <atomic>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <limits>
#include <zlib.h>
#include <tinyxml2.h>
#include <omp.h>
#include <iomanip>
#include <map>

namespace fs = std::filesystem;
using namespace tinyxml2;

namespace {
bool verbose() {
    static const bool v = []() {
        const char* e = std::getenv("CUDA_DEMUX_VERBOSE");
        return e && e[0] && e[0] != '0';
    }();
    return v;
}
#define VLOG if (!verbose()) {} else std::cout
}

// Represents the structure of a sequencing run
struct RunStructure {
    int read1_cycles = 0;
    int index1_cycles = 0;
    int index2_cycles = 0;
    int read2_cycles = 0;
    std::vector<int> read_segments; // 0:R1, 1:I1, 2:I2, 3:R2
};

// Forward declarations for parsing functions
std::vector<LaneBclData> parse_legacy_bcl(const fs::path& bcl_dir, const RunStructure& run_structure, int total_cycles);
std::vector<LaneBclData> parse_cbcl(const fs::path& bcl_dir, const RunStructure& run_structure, int total_cycles);

static void apply_run_structure(LaneBclData& lane, const RunStructure& rs, int total_cycles) {
    lane.total_cycles = total_cycles;
    lane.r1_len = rs.read1_cycles;
    lane.i1_len = rs.index1_cycles;
    lane.i2_len = rs.index2_cycles;
    lane.r2_len = rs.read2_cycles;
    lane.read_segments = rs.read_segments;
}

// Main C++ function to parse BCL data
std::vector<LaneBclData> parse_bcl(const std::string& bcl_folder) {
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
    fs::path basecalls_dir = bcl_dir / "Data" / "Intensities" / "BaseCalls";
    
    if (!fs::exists(basecalls_dir)) {
        std::cerr << "Error: BaseCalls directory not found: " << basecalls_dir.string() << std::endl;
        return {};
    }
    
    VLOG << "Scanning directory: " << basecalls_dir.string() << std::endl;
    
    bool has_cbcl = false;
    bool has_legacy_bcl = false;
    std::vector<fs::path> cbcl_files;
    std::vector<fs::path> legacy_bcl_files;
    
    // First, scan the BaseCalls directory (may contain files in some runs)
    for(const auto& entry : fs::directory_iterator(basecalls_dir)) {
        VLOG << "  Found: " << entry.path().filename().string() << std::endl;
        if (entry.path().extension() == ".cbcl") {
            has_cbcl = true;
            cbcl_files.push_back(entry.path());
        } else if (entry.path().extension() == ".gz" && entry.path().stem().extension() == ".bcl") {
            has_legacy_bcl = true;
            legacy_bcl_files.push_back(entry.path());
        }
    }

    // Then, scan lane directories for cycle subfolders and files
    VLOG << "Scanning lane directories for BCL/CBCL files..." << std::endl;
    for(const auto& lane_entry : fs::directory_iterator(basecalls_dir)) {
        const std::string lane_name = lane_entry.path().filename().string();
        if (lane_entry.is_directory() && lane_name.rfind("L00", 0) == 0) {
            VLOG << "  Scanning lane: " << lane_entry.path().filename().string() << std::endl;
            // Direct files in lane dir
            for (const auto& entry : fs::directory_iterator(lane_entry.path())) {
                if (entry.path().extension() == ".cbcl") { has_cbcl = true; }
                else if (entry.path().extension() == ".gz" && entry.path().stem().extension() == ".bcl") { has_legacy_bcl = true; }
            }
            // Cycle directories
            for(const auto& entry : fs::directory_iterator(lane_entry.path())) {
                const std::string cyc_name = entry.path().filename().string();
                if (entry.is_directory() && (!cyc_name.empty() && cyc_name[0] == 'C')) {
                    for(const auto& subentry : fs::directory_iterator(entry.path())) {
                        if (subentry.path().extension() == ".cbcl") { has_cbcl = true; }
                        else if (subentry.path().extension() == ".gz" && subentry.path().stem().extension() == ".bcl") { has_legacy_bcl = true; }
                    }
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

std::vector<LaneBclData> parse_legacy_bcl(const fs::path& bcl_dir, const RunStructure& run_structure, int total_cycles) {
    fs::path basecalls_dir = bcl_dir / "Data" / "Intensities" / "BaseCalls";
    std::cout << "Starting legacy BCL parsing in directory: " << basecalls_dir.string() << std::endl;

    std::vector<fs::path> lane_dirs;
    for (const auto& entry : fs::directory_iterator(basecalls_dir)) {
        if (fs::is_directory(entry) && entry.path().filename().string().find("L00") == 0) {
            lane_dirs.push_back(entry.path());
        }
    }
    if (lane_dirs.empty()) {
        lane_dirs.push_back(basecalls_dir);
    }
    std::sort(lane_dirs.begin(), lane_dirs.end());
    std::cout << "Found " << lane_dirs.size() << " lanes to process (legacy BCL)" << std::endl;

    std::vector<LaneBclData> all_lanes;

    for (const auto& lane_dir : lane_dirs) {
        std::string lane_name = lane_dir.filename().string();
        std::cout << "Processing lane: " << lane_name << std::endl;

        LaneBclData lane;
        apply_run_structure(lane, run_structure, total_cycles);
        try {
            if (lane_name.size() >= 4 && lane_name[0] == 'L') {
                lane.lane = std::stoi(lane_name.substr(1));
            }
        } catch (...) { lane.lane = 1; }

        lane.bcl.assign(total_cycles, {});
        size_t num_clusters = 0;

        for (int c = 1; c <= total_cycles; ++c) {
            fs::path cycle_dir = lane_dir / ("C" + std::to_string(c) + ".1");
            fs::path bcl_file = cycle_dir / (lane_name + "_1.bcl.gz");
            if (!fs::exists(bcl_file)) {
                bcl_file = cycle_dir / "s_1_1101.bcl.gz";
                if (!fs::exists(bcl_file)) {
                    throw std::runtime_error("BCL file not found for cycle " + std::to_string(c) +
                                             " in " + cycle_dir.string());
                }
            }
            uint32_t current_clusters = 0;
            std::vector<char> data = read_bcl_gz_file(bcl_file, current_clusters);
            if (c == 1) {
                num_clusters = current_clusters;
            } else if (current_clusters != num_clusters) {
                throw std::runtime_error("Inconsistent cluster count across BCL files for lane " +
                                         lane_name);
            }
            lane.bcl[c - 1].assign(data.begin(), data.end());
        }

        lane.num_clusters = num_clusters;
        std::cout << "Loaded " << total_cycles << " BCL files for " << num_clusters
                  << " clusters in " << lane_name << "." << std::endl;
        all_lanes.push_back(std::move(lane));
    }
    return all_lanes;
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

// CBCL parsing functions based on actual hex dump analysis
// Diagnostics and alternate decode modes for platform variance
static inline char base_from2(uint8_t b2) {
    static const char MAP[4] = {'A','C','G','T'};
    return MAP[b2 & 0x3];
}

struct DecodeStats {
    size_t total = 0;
    size_t valid = 0;
    size_t acgt[4] = {0,0,0,0};
    size_t q_nonzero = 0;
    size_t q_zero = 0;
};

static DecodeStats tally_decode(const std::vector<uint8_t>& bases_idx,
                                const std::vector<uint8_t>& quals) {
    DecodeStats s; s.total = bases_idx.size(); s.valid = bases_idx.size();
    for (size_t i=0; i<bases_idx.size(); ++i) {
        char b = base_from2(bases_idx[i]);
        if (b=='A') s.acgt[0]++; else if (b=='C') s.acgt[1]++; else if (b=='G') s.acgt[2]++; else if (b=='T') s.acgt[3]++;
        if (i < quals.size()) { if (quals[i] > 0) s.q_nonzero++; else s.q_zero++; }
    }
    return s;
}

static void print_decode_stats(const char* label, const DecodeStats& s) {
    size_t acgt_total = s.acgt[0]+s.acgt[1]+s.acgt[2]+s.acgt[3];
    double pct = s.valid ? (100.0 * (double)acgt_total / (double)s.valid) : 0.0;
    std::cout << "    " << label << ": total=" << s.total
              << " A:" << s.acgt[0] << " C:" << s.acgt[1]
              << " G:" << s.acgt[2] << " T:" << s.acgt[3]
              << " q>0:" << s.q_nonzero << " q=0:" << s.q_zero
              << "  ACGT%=" << pct << std::endl;
}

static void decode_mode_per_cluster_byte(const std::vector<uint8_t>& uncmp, uint32_t n,
                                         std::vector<uint8_t>& out_bases, std::vector<uint8_t>& out_quals) {
    out_bases.clear(); out_quals.clear();
    uint32_t limit = std::min<uint32_t>(n, (uint32_t)uncmp.size());
    out_bases.reserve(limit); out_quals.reserve(limit);
    for (uint32_t i=0; i<limit; ++i) {
        uint8_t b = uncmp[i];
        out_bases.push_back(b & 0x03);
        out_quals.push_back((b >> 2) & 0x3F);
    }
}

static void decode_mode_interleaved_nibbles(const std::vector<uint8_t>& uncmp, uint32_t n,
                                            std::vector<uint8_t>& out_bases, std::vector<uint8_t>& out_quals) {
    out_bases.clear(); out_quals.clear(); out_bases.reserve(n); out_quals.reserve(n);
    for (uint32_t i=0; i<n; ++i) {
        uint32_t byte_index = i >> 1; if (byte_index >= uncmp.size()) break;
        uint8_t byte = uncmp[byte_index];
        uint8_t nib = (i & 1) ? ((byte >> 4) & 0x0F) : (byte & 0x0F);
        uint8_t bi = nib & 0x03; uint8_t q2 = (nib >> 2) & 0x03;
        out_bases.push_back(bi);
        out_quals.push_back((uint8_t)(q2*10+2));
    }
}

static void decode_mode_separate_streams(const std::vector<uint8_t>& uncmp, uint32_t n,
                                         int bits_per_base, int bits_per_q,
                                         std::vector<uint8_t>& out_bases, std::vector<uint8_t>& out_quals) {
    out_bases.assign(n, 0); out_quals.assign(n, 0);
    const size_t baseBits = (size_t)bits_per_base * n; // expect 2*n
    const size_t baseBytes = (baseBits + 7) >> 3;
    if (uncmp.size() < baseBytes) return;
    size_t bitpos = 0;
    for (uint32_t i=0; i<n; ++i) {
        size_t byteIdx = bitpos >> 3; int shift = (int)(bitpos & 7);
        uint8_t b = uncmp[byteIdx]; uint8_t val;
        if (shift <= 6) val = (b >> shift) & 0x03;
        else {
            if (byteIdx + 1 >= uncmp.size()) { val = 0; }
            else { uint16_t w = (uint16_t)b | ((uint16_t)uncmp[byteIdx+1] << 8); val = (w >> shift) & 0x03; }
        }
        out_bases[i] = val; bitpos += bits_per_base;
    }
    const uint8_t* qptr = (baseBytes < uncmp.size()) ? &uncmp[baseBytes] : nullptr;
    if (!qptr) return;
    int qBytesPer = (bits_per_q + 7) >> 3; if (qBytesPer <= 0) qBytesPer = 1;
    size_t need = (size_t)qBytesPer * n; if (uncmp.size() < baseBytes + need) return;
    for (uint32_t i=0; i<n; ++i) {
        uint8_t qraw = qptr[i * qBytesPer];
        out_quals[i] = (bits_per_q >= 6) ? (qraw & 0x3F) : qraw;
    }
}
// Maintain minimal global decode format inferred from the header of the current CBCL file.
static int g_cbcl_bits_per_base = 2;
static int g_cbcl_bits_per_q = 2;
CbclHeader parse_cbcl_header(std::ifstream& file) {
    CbclHeader header;
    
    // Based on the hex dump and NovaSeqX CBCL format:
    // Bytes 0-1: Version (01 00 = version 1)
    // Bytes 2-5: Header size (31 06 00 00 = 1585 bytes)
    // Byte 6: Bits per basecall (02 = 2 bits)
    // Byte 7: Bits per q-score (02 = 2 bits)
    // Bytes 8-11: Number of quality bins (04 00 00 00 = 4)
    // Bytes 12-15: Quality bin remapping (00 00 00 00)
    // Bytes 16-19: Unknown field (00 00 00 00)
    // Bytes 20-23: Unknown field (01 00 00 00)
    // Bytes 24-27: Unknown field (0C 00 00 00)
    // Bytes 28-31: Unknown field (02 00 00 00)
    // Bytes 32-35: Unknown field (18 00 00 00)
    // Bytes 36-39: Unknown field (03 00 00 00)
    // Bytes 40-43: Unknown field (28 00 00 00)
    // Bytes 44-47: Number of tiles (60 00 00 00 = 96)
    
    // Read version as 2 bytes
    uint16_t version;
    file.read(reinterpret_cast<char*>(&version), sizeof(uint16_t));
    header.version = version;
    
    // Read header size
    file.read(reinterpret_cast<char*>(&header.header_size), sizeof(uint32_t));
    
    // Read bits per basecall and quality
    uint8_t bits_per_basecall, bits_per_quality;
    file.read(reinterpret_cast<char*>(&bits_per_basecall), 1);
    file.read(reinterpret_cast<char*>(&bits_per_quality), 1);
    header.bits_per_basecall = bits_per_basecall;
    header.bits_per_quality = bits_per_quality;
    // Expose to parser state for subsequent block decoding
    g_cbcl_bits_per_base = static_cast<int>(bits_per_basecall);
    g_cbcl_bits_per_q = static_cast<int>(bits_per_quality);
    
    // Read number of quality bins (bins for binning quality scores)
    file.read(reinterpret_cast<char*>(&header.num_quality_bins), sizeof(uint32_t));
    
    // Skip quality score remapping table (num_quality_bins entries)
    file.seekg(header.num_quality_bins, std::ios::cur);
    
    // Skip unknown fields to get to tile count at offset 0x2C (44)
    // Current position is at byte 16, need to get to byte 44
    file.seekg(28, std::ios::cur);  // Skip 28 bytes
    
    // Read number of tiles
    file.read(reinterpret_cast<char*>(&header.num_tiles), sizeof(uint32_t));
    
    // Set other fields based on CBCL format
    header.compression_type = 2;  // Gzip compression
    header.num_cycles = 1;  // CBCL files store one cycle at a time
    header.bits_per_filter = 1;  // Standard for filter values
    header.num_bases_per_cycle = 1;  // One base per cycle
    
    // Calculate tile list offset - it starts right after the current position
    header.tile_list_offset = file.tellg();
    
    VLOG << "CBCL Header: version=" << header.version
         << ", header_size=" << header.header_size
         << ", bits_per_basecall=" << (int)header.bits_per_basecall
         << ", bits_per_quality=" << (int)header.bits_per_quality
         << ", num_quality_bins=" << header.num_quality_bins
         << ", num_tiles=" << header.num_tiles
         << ", tile_list_offset=0x" << std::hex << header.tile_list_offset << std::dec << std::endl;
    
    return header;
}

std::vector<CbclTileInfo> parse_cbcl_tile_list(std::ifstream& file, const CbclHeader& header) {
    std::vector<CbclTileInfo> tiles;
    
    file.seekg(header.tile_list_offset);
    
    VLOG << "Reading " << header.num_tiles << " tile records from offset 0x"
         << std::hex << header.tile_list_offset << std::dec << std::endl;
    
    // Based on CBCL format specification, each tile record contains:
    // - Tile number (4 bytes)
    // - Number of clusters (4 bytes)
    // - Uncompressed block size (4 bytes)
    // - Compressed block size (4 bytes)
    
    for (uint32_t i = 0; i < header.num_tiles; ++i) {
        CbclTileInfo tile;
        
        // Read tile number
        file.read(reinterpret_cast<char*>(&tile.tile_id), sizeof(uint32_t));
        
        // Read number of clusters
        file.read(reinterpret_cast<char*>(&tile.num_clusters), sizeof(uint32_t));
        
        // Read uncompressed block size
        file.read(reinterpret_cast<char*>(&tile.uncompressed_block_size), sizeof(uint32_t));
        
        // Read compressed block size
        file.read(reinterpret_cast<char*>(&tile.compressed_block_size), sizeof(uint32_t));
        
        // File offset will be calculated cumulatively after we read all tiles
        tile.file_offset = 0;  // Will be set later
        
        tiles.push_back(tile);
        
        VLOG << "    Tile " << std::dec << tile.tile_id
             << ": clusters=" << tile.num_clusters
             << ", uncompressed=" << tile.uncompressed_block_size
             << ", compressed=" << tile.compressed_block_size << std::endl;
    }

    uint64_t current_offset = header.header_size;
    for (auto& tile : tiles) {
        tile.file_offset = current_offset;
        current_offset += tile.compressed_block_size;
    }

    VLOG << "Parsed " << tiles.size() << " tile entries" << std::endl;
    return tiles;
}

std::vector<uint8_t> decompress_cbcl_data(const std::vector<char>& compressed_data, uint32_t uncompressed_size) {
    std::vector<uint8_t> uncompressed_data(uncompressed_size);
    
    // Try different zlib decompression methods
    z_stream strm = {};
    
    // Method 1: Try raw deflate (no header)
    strm.avail_in = compressed_data.size();
    strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(compressed_data.data()));
    
    int result = inflateInit2(&strm, MAX_WBITS); // Raw deflate
    if (result == Z_OK) {
        uLong actual_size = uncompressed_size;
        result = uncompress(reinterpret_cast<Bytef*>(uncompressed_data.data()), &actual_size,
                           reinterpret_cast<const Bytef*>(compressed_data.data()), compressed_data.size());
        
        if (result == Z_OK) {
            uncompressed_data.resize(actual_size);
            inflateEnd(&strm);
            return uncompressed_data;
        }
        inflateEnd(&strm);
    }
    
    // Method 2: Try gzip format
    strm = {};
    strm.avail_in = compressed_data.size();
    strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(compressed_data.data()));
    
    result = inflateInit2(&strm, 16 + MAX_WBITS); // Gzip format
    if (result == Z_OK) {
        uLong actual_size = uncompressed_size;
        result = uncompress(reinterpret_cast<Bytef*>(uncompressed_data.data()), &actual_size,
                           reinterpret_cast<const Bytef*>(compressed_data.data()), compressed_data.size());
        
        if (result == Z_OK) {
            uncompressed_data.resize(actual_size);
            inflateEnd(&strm);
            return uncompressed_data;
        }
        inflateEnd(&strm);
    }
    
    // Method 3: Try streaming decompression with larger buffer
    strm = {};
    strm.avail_in = compressed_data.size();
    strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(compressed_data.data()));
    
    result = inflateInit2(&strm, MAX_WBITS); // Raw deflate
    if (result == Z_OK) {
        std::vector<uint8_t> temp_buffer(uncompressed_size * 2); // Larger buffer
        strm.avail_out = temp_buffer.size();
        strm.next_out = temp_buffer.data();
        
        result = inflate(&strm, Z_FINISH);
        if (result == Z_STREAM_END) {
            uncompressed_data.resize(strm.total_out);
            std::copy(temp_buffer.begin(), temp_buffer.begin() + strm.total_out, uncompressed_data.begin());
            inflateEnd(&strm);
            return uncompressed_data;
        }
        inflateEnd(&strm);
    }
    
    // Method 4: Try without specifying uncompressed size (let zlib determine it)
    strm = {};
    strm.avail_in = compressed_data.size();
    strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(compressed_data.data()));
    
    result = inflateInit2(&strm, MAX_WBITS);
    if (result == Z_OK) {
        std::vector<uint8_t> temp_buffer(32768); // Start with reasonable size
        int ret;
        do {
            strm.avail_out = temp_buffer.size() - strm.total_out;
            strm.next_out = temp_buffer.data() + strm.total_out;
            ret = inflate(&strm, Z_NO_FLUSH);
            
            if (ret != Z_OK && ret != Z_STREAM_END) {
                break;
            }
            
            if (strm.avail_out == 0 && ret != Z_STREAM_END) {
                temp_buffer.resize(temp_buffer.size() * 2);
            }
        } while (ret != Z_STREAM_END);
        
        if (ret == Z_STREAM_END) {
            uncompressed_data.resize(strm.total_out);
            std::copy(temp_buffer.begin(), temp_buffer.begin() + strm.total_out, uncompressed_data.begin());
            inflateEnd(&strm);
            return uncompressed_data;
        }
        inflateEnd(&strm);
    }
    
    std::cerr << "Error: Failed to decompress CBCL data. zlib error: " << result 
              << " (Z_DATA_ERROR = -3, Z_BUF_ERROR = -5)" << std::endl;
    return {};
}

CbclBlock parse_cbcl_block(std::ifstream& file, const CbclTileInfo& tile_info) {
    CbclBlock block;
    
    // Seek to the tile data
    file.seekg(tile_info.file_offset);
    
    // Read the compressed data
    std::vector<char> compressed_data(tile_info.compressed_block_size);
    file.read(compressed_data.data(), tile_info.compressed_block_size);
    
    if (file.gcount() != static_cast<std::streamsize>(tile_info.compressed_block_size)) {
        std::cerr << "      Error: Could only read " << file.gcount() 
                  << " bytes, expected " << tile_info.compressed_block_size << std::endl;
        return block;
    }
    
    // Decompress the data using zlib
    std::vector<uint8_t> uncompressed_data(tile_info.uncompressed_block_size);
    
    z_stream strm = {};
    strm.avail_in = tile_info.compressed_block_size;
    strm.next_in = reinterpret_cast<Bytef*>(compressed_data.data());
    strm.avail_out = tile_info.uncompressed_block_size;
    strm.next_out = uncompressed_data.data();
    
    // Initialize for gzip decompression (16 + MAX_WBITS for gzip format)
    int ret = inflateInit2(&strm, 16 + MAX_WBITS);
    if (ret != Z_OK) {
        std::cerr << "      Error: Failed to initialize zlib for decompression" << std::endl;
        return block;
    }
    
    ret = inflate(&strm, Z_FINISH);
    if (ret != Z_STREAM_END) {
        std::cerr << "      Error: Failed to decompress data, zlib error: " << ret << std::endl;
        inflateEnd(&strm);
        return block;
    }
    
    uint32_t decompressed_size = strm.total_out;
    inflateEnd(&strm);
    
    VLOG << "      Decompressed " << decompressed_size << " bytes from "
         << tile_info.compressed_block_size << " compressed bytes" << std::endl;
    
    // Parse the decompressed data into basecalls + qualities per cluster.
    std::vector<uint8_t> basecalls; basecalls.reserve(tile_info.num_clusters);
    std::vector<uint8_t> qualities; qualities.reserve(tile_info.num_clusters);
    std::vector<uint8_t> filters;

    bool do_test = (std::getenv("CUDA_DEMUX_CBCL_TEST") != nullptr);
    if (do_test) {
        std::vector<uint8_t> b_pc, q_pc, b_nib, q_nib, b_sep, q_sep;
        decode_mode_per_cluster_byte(uncompressed_data, tile_info.num_clusters, b_pc, q_pc);
        decode_mode_interleaved_nibbles(uncompressed_data, tile_info.num_clusters, b_nib, q_nib);
        decode_mode_separate_streams(uncompressed_data, tile_info.num_clusters, 2, g_cbcl_bits_per_q, b_sep, q_sep);
        auto sc = [](const DecodeStats& s){ return (long long)(s.acgt[0]+s.acgt[1]+s.acgt[2]+s.acgt[3]); };
        DecodeStats spc = tally_decode(b_pc, q_pc);
        DecodeStats snib = tally_decode(b_nib, q_nib);
        DecodeStats ssep = tally_decode(b_sep, q_sep);
        std::cout << "  [CBCL TEST] tile " << tile_info.tile_id << ": clusters=" << tile_info.num_clusters << std::endl;
        print_decode_stats("per-cluster-byte", spc);
        print_decode_stats("nibble-interlv   ", snib);
        print_decode_stats("separate-stream  ", ssep);
        long long cpc=sc(spc), cnib=sc(snib), csep=sc(ssep);
        if (csep >= cpc && csep >= cnib && !b_sep.empty()) { basecalls=std::move(b_sep); qualities=std::move(q_sep); }
        else if (cpc >= cnib && !b_pc.empty()) { basecalls=std::move(b_pc); qualities=std::move(q_pc); }
        else { basecalls=std::move(b_nib); qualities=std::move(q_nib); }
    } else {
        // Heuristic: if qualities are bytes (>=6 bits), prefer separate streams first; else nibble fallback
        if (g_cbcl_bits_per_q >= 6) {
            decode_mode_separate_streams(uncompressed_data, tile_info.num_clusters, 2, g_cbcl_bits_per_q, basecalls, qualities);
            if (basecalls.empty()) {
                decode_mode_per_cluster_byte(uncompressed_data, tile_info.num_clusters, basecalls, qualities);
            }
        }
        if (basecalls.empty()) {
            decode_mode_interleaved_nibbles(uncompressed_data, tile_info.num_clusters, basecalls, qualities);
        }
    }
    
    VLOG << "      Extracted " << basecalls.size() << " basecalls and "
         << qualities.size() << " quality scores" << std::endl;
    
    // Store the extracted data
    block.num_clusters = tile_info.num_clusters;
    block.basecalls = basecalls;
    block.qualities = qualities;
    block.filters = filters; // Empty for now, as filter data is in separate files
    
    return block;
}

std::vector<LaneBclData> parse_cbcl(const fs::path& bcl_dir, const RunStructure& run_structure, int total_cycles) {
    fs::path basecalls_dir = bcl_dir / "Data" / "Intensities" / "BaseCalls";

    std::vector<fs::path> lane_dirs;
    for (const auto& entry : fs::directory_iterator(basecalls_dir)) {
        if (fs::is_directory(entry) && entry.path().filename().string().find("L00") == 0) {
            lane_dirs.push_back(entry.path());
        }
    }

    if (lane_dirs.empty()) {
        lane_dirs.push_back(basecalls_dir);
    }

    std::sort(lane_dirs.begin(), lane_dirs.end());
    std::cout << "Found " << lane_dirs.size() << " lanes to process" << std::endl;

    std::vector<LaneBclData> all_lanes;

    for (const auto& lane_dir : lane_dirs) {
        std::cout << "Processing lane: " << lane_dir.string() << std::endl;

        int lane_number = 1;
        try {
            std::string lname = lane_dir.filename().string();
            if (lname.size() >= 4 && lname[0] == 'L') {
                lane_number = std::stoi(lname.substr(1));
            }
        } catch (...) { lane_number = 1; }

        // Per-tile filter files are named s_<lane>_<tile>.filter on every Illumina platform.
        const std::string filter_prefix = "s_" + std::to_string(lane_number) + "_";
        std::vector<fs::path> filter_files;

        for (const auto& entry : fs::directory_iterator(lane_dir)) {
            const std::string fname = entry.path().filename().string();
            if (entry.path().extension() == ".filter" &&
                fname.rfind(filter_prefix, 0) == 0) {
                filter_files.push_back(entry.path());
            }
        }

        if (filter_files.empty()) {
            fs::path single_filter_file = lane_dir / ("s_" + std::to_string(lane_number) + ".filter");
            if (fs::exists(single_filter_file)) {
                filter_files.push_back(single_filter_file);
            }
        }

        if (filter_files.empty()) {
            std::cerr << "Error: No filter files found in " << lane_dir.string() << std::endl;
            std::cerr << "Available files in directory:" << std::endl;
            for(const auto& entry : fs::directory_iterator(lane_dir)) {
                std::cerr << "  - " << entry.path().filename().string() << std::endl;
            }
            throw std::runtime_error("No filter files found in " + lane_dir.string());
        }

        std::cout << "Found " << filter_files.size() << " filter files" << std::endl;
        std::sort(filter_files.begin(), filter_files.end());

        std::vector<uint8_t> passing_clusters;
        uint32_t num_clusters_total = 0;
        uint32_t num_clusters_passed = 0;

        for (const auto& filter_file : filter_files) {
            std::ifstream filter_stream(filter_file, std::ios::binary);
            if (!filter_stream) {
                throw std::runtime_error("Failed to open filter file: " + filter_file.string());
            }

            uint32_t filter_version, unknown_field, tile_clusters;
            filter_stream.read(reinterpret_cast<char*>(&filter_version), sizeof(filter_version));

            if (filter_version == 0) {
                filter_stream.read(reinterpret_cast<char*>(&unknown_field), sizeof(unknown_field));
                filter_stream.read(reinterpret_cast<char*>(&tile_clusters), sizeof(tile_clusters));
            } else {
                filter_stream.read(reinterpret_cast<char*>(&tile_clusters), sizeof(tile_clusters));
            }

            VLOG << "  Reading " << filter_file.filename().string()
                 << ": version=" << filter_version << ", clusters=" << tile_clusters << std::endl;

            const size_t base = passing_clusters.size();
            passing_clusters.resize(base + tile_clusters);
            filter_stream.read(reinterpret_cast<char*>(passing_clusters.data() + base), tile_clusters);
            if (filter_stream.gcount() != static_cast<std::streamsize>(tile_clusters)) {
                throw std::runtime_error("Truncated filter file: " + filter_file.string());
            }
            uint32_t tile_passed = 0;
            for (size_t i = 0; i < tile_clusters; ++i) {
                if (passing_clusters[base + i] == 1) {
                    passing_clusters[base + i] = 1;
                    ++tile_passed;
                } else {
                    passing_clusters[base + i] = 0;
                }
            }
            num_clusters_passed += tile_passed;
            num_clusters_total += tile_clusters;

            VLOG << "    Tile " << filter_file.filename().string()
                 << ": " << tile_passed << "/" << tile_clusters << " clusters passed QC" << std::endl;
        }

        std::cout << "Combined filter data: " << num_clusters_total << " total clusters, "
                  << num_clusters_passed << " passing QC." << std::endl;

        std::vector<fs::path> cbcl_files;
        for (const auto& entry : fs::directory_iterator(lane_dir)) {
            if (entry.path().extension() == ".cbcl") {
                cbcl_files.push_back(entry.path());
            }
        }
        for (const auto& entry : fs::directory_iterator(lane_dir)) {
            const std::string cyc_name = entry.path().filename().string();
            if (entry.is_directory() && !cyc_name.empty() && cyc_name[0] == 'C') {
                for (const auto& subentry : fs::directory_iterator(entry.path())) {
                    if (subentry.path().extension() == ".cbcl") {
                        cbcl_files.push_back(subentry.path());
                    }
                }
            }
        }
        std::sort(cbcl_files.begin(), cbcl_files.end());
        std::cout << "Found " << cbcl_files.size() << " CBCL files" << std::endl;
        if (cbcl_files.empty()) {
            throw std::runtime_error("No CBCL files found in " + lane_dir.string());
        }

        if (num_clusters_passed == 0) {
            std::cerr << "Warning: No clusters passed QC filter in lane "
                      << lane_dir.filename().string() << std::endl;
            continue;
        }
        std::cout << "After QC filtering, keeping " << num_clusters_passed << " clusters ("
                  << (num_clusters_total - num_clusters_passed) << " filtered out)" << std::endl;

        std::map<int, std::vector<fs::path>> cycle_to_cbcl_files;
        for (const auto& cbcl_path : cbcl_files) {
            std::string parent_dir = cbcl_path.parent_path().filename().string();
            if (!parent_dir.empty() && parent_dir[0] == 'C' &&
                parent_dir.find('.') != std::string::npos) {
                int cycle = std::stoi(parent_dir.substr(1, parent_dir.find('.') - 1));
                cycle_to_cbcl_files[cycle].push_back(cbcl_path);
            }
        }
        for (auto& [cycle, files] : cycle_to_cbcl_files) {
            std::sort(files.begin(), files.end());
        }
        std::cout << "Found CBCL files for " << cycle_to_cbcl_files.size() << " cycles" << std::endl;

        // pass-index map (uint32_t halves the memory vs int)
        constexpr uint32_t kFiltered = std::numeric_limits<uint32_t>::max();
        std::vector<uint32_t> pass_index_map(num_clusters_total, kFiltered);
        {
            uint32_t compact_idx = 0;
            for (uint32_t i = 0; i < num_clusters_total; ++i) {
                if (passing_clusters[i]) {
                    pass_index_map[i] = compact_idx++;
                }
            }
        }
        passing_clusters.clear();
        passing_clusters.shrink_to_fit();

        LaneBclData lane;
        apply_run_structure(lane, run_structure, total_cycles);
        lane.lane = lane_number;
        lane.num_clusters = num_clusters_passed;
        lane.bcl.assign(total_cycles, std::vector<uint8_t>(num_clusters_passed));

        std::cout << "Reading " << cycle_to_cbcl_files.size() << " cycles using up to "
                  << omp_get_max_threads() << " threads for decompression..." << std::endl;

        std::vector<int> cycles_to_process;
        cycles_to_process.reserve(cycle_to_cbcl_files.size());
        for (auto& [cyc, _] : cycle_to_cbcl_files) {
            cycles_to_process.push_back(cyc);
        }

        std::atomic<bool> hit_error{false};
        std::string first_error;
        #pragma omp critical
        first_error.clear();

        #pragma omp parallel for schedule(dynamic)
        for (int ci = 0; ci < static_cast<int>(cycles_to_process.size()); ++ci) {
            if (hit_error.load()) continue;
            int cycle = cycles_to_process[ci];
            const auto& cbcl_paths = cycle_to_cbcl_files[cycle];
            try {
                uint32_t cluster_offset = 0;
                for (const auto& cbcl_path : cbcl_paths) {
                    std::ifstream cbcl_file(cbcl_path, std::ios::binary);
                    if (!cbcl_file) {
                        throw std::runtime_error("Failed to open CBCL file: " + cbcl_path.string());
                    }
                    CbclHeader header = parse_cbcl_header(cbcl_file);
                    std::vector<CbclTileInfo> tiles = parse_cbcl_tile_list(cbcl_file, header);

                    for (const auto& tile : tiles) {
                        CbclBlock block = parse_cbcl_block(cbcl_file, tile);
                        if (block.basecalls.empty()) continue;
                        const int cycle_index = cycle - 1;
                        const uint32_t n = static_cast<uint32_t>(block.basecalls.size());
                        for (uint32_t i = 0; i < n; ++i) {
                            uint32_t global_idx = cluster_offset + i;
                            if (global_idx >= pass_index_map.size()) break;
                            uint32_t out_idx = pass_index_map[global_idx];
                            if (out_idx == kFiltered) continue;
                            uint8_t basecall = block.basecalls[i] & 0x03;
                            uint8_t quality = (i < block.qualities.size()) ? block.qualities[i] : 0;
                            lane.bcl[cycle_index][out_idx] =
                                static_cast<uint8_t>((quality << 2) | basecall);
                        }
                        cluster_offset += tile.num_clusters;
                    }
                }
            } catch (const std::exception& e) {
                bool expected = false;
                if (hit_error.compare_exchange_strong(expected, true)) {
                    #pragma omp critical
                    first_error = std::string("cycle ") + std::to_string(cycle) +
                                  ": " + e.what();
                }
            }
        }

        if (hit_error.load()) {
            throw std::runtime_error("CBCL parse failed: " + first_error);
        }

        std::cout << "Lane " << lane_dir.filename().string()
                  << ": prepared " << num_clusters_passed << " clusters across "
                  << total_cycles << " cycles." << std::endl;
        all_lanes.push_back(std::move(lane));
    }

    std::cout << "Total lanes prepared: " << all_lanes.size() << std::endl;
    return all_lanes;
}
