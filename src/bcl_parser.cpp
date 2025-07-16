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
#include <iomanip>
#include <map>

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

// CBCL parsing functions based on actual hex dump analysis
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
    
    std::cout << "CBCL Header: version=" << header.version 
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
    
    std::cout << "Reading " << header.num_tiles << " tile records from offset 0x" 
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
        
        std::cout << "    Tile " << std::dec << tile.tile_id 
                  << ": clusters=" << tile.num_clusters
                  << ", uncompressed=" << tile.uncompressed_block_size
                  << ", compressed=" << tile.compressed_block_size << std::endl;
    }
    
    // Now calculate the actual file offsets for each tile's compressed data
    // The data starts immediately after the header
    uint64_t current_offset = header.header_size;
    
    for (auto& tile : tiles) {
        tile.file_offset = current_offset;
        current_offset += tile.compressed_block_size;
        
        std::cout << "    Tile " << tile.tile_id 
                  << " data at offset 0x" << std::hex << tile.file_offset << std::dec << std::endl;
    }
    
    std::cout << "Parsed " << tiles.size() << " tile entries" << std::endl;
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
    
    std::cout << "      Decompressed " << decompressed_size << " bytes from " 
              << tile_info.compressed_block_size << " compressed bytes" << std::endl;
    
    // Parse the decompressed data
    // CBCL format packs data as:
    // 1. Basecalls (2 bits per base, packed 4 per byte)
    // 2. Quality scores (variable bits per quality, remapped via quality bins)
    // 3. Filter data (if present)
    
    uint32_t offset = 0;
    std::vector<uint8_t> basecalls;
    std::vector<uint8_t> qualities;
    std::vector<uint8_t> filters;
    
    // Extract basecalls (2 bits per basecall)
    uint32_t basecall_bytes = (tile_info.num_clusters + 3) / 4; // Round up
    basecalls.reserve(tile_info.num_clusters);
    
    for (uint32_t i = 0; i < tile_info.num_clusters && offset < decompressed_size; ++i) {
        uint32_t byte_index = offset + (i / 4);
        uint32_t bit_offset = (i % 4) * 2;
        
        if (byte_index < decompressed_size) {
            uint8_t byte = uncompressed_data[byte_index];
            uint8_t basecall = (byte >> bit_offset) & 0x03;
            basecalls.push_back(basecall);
        }
    }
    offset += basecall_bytes;
    
    // Extract quality scores
    // For NovaSeqX with 2 bits per quality, we have 4 quality values per byte
    uint32_t quality_bytes = (tile_info.num_clusters + 3) / 4; // 2 bits per quality
    qualities.reserve(tile_info.num_clusters);
    
    for (uint32_t i = 0; i < tile_info.num_clusters && offset < decompressed_size; ++i) {
        uint32_t byte_index = offset + (i / 4);
        uint32_t bit_offset = (i % 4) * 2;
        
        if (byte_index < decompressed_size) {
            uint8_t byte = uncompressed_data[byte_index];
            uint8_t quality = (byte >> bit_offset) & 0x03;
            // Map 2-bit quality to standard quality scores
            // This is a simplified mapping - real implementation would use quality bins from header
            uint8_t mapped_quality = quality * 10 + 2; // Maps 0-3 to 2, 12, 22, 32
            qualities.push_back(mapped_quality);
        }
    }
    
    std::cout << "      Extracted " << basecalls.size() << " basecalls and " 
              << qualities.size() << " quality scores" << std::endl;
    
    // Store the extracted data
    block.num_clusters = tile_info.num_clusters;
    block.basecalls = basecalls;
    block.qualities = qualities;
    block.filters = filters; // Empty for now, as filter data is in separate files
    
    return block;
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
        
        uint32_t filter_version, unknown_field, tile_clusters;
        filter_stream.read(reinterpret_cast<char*>(&filter_version), sizeof(filter_version));
        
        // NovaSeqX filter format has an additional field after version
        if (filter_version == 0) {
            // Read unknown field (appears to be 3 in test data)
            filter_stream.read(reinterpret_cast<char*>(&unknown_field), sizeof(unknown_field));
            // Read actual cluster count
            filter_stream.read(reinterpret_cast<char*>(&tile_clusters), sizeof(tile_clusters));
        } else {
            // Standard filter format - cluster count follows version
            filter_stream.read(reinterpret_cast<char*>(&tile_clusters), sizeof(tile_clusters));
        }
        
        std::cout << "  Reading " << filter_file.filename().string() 
                  << ": version=" << filter_version << ", clusters=" << tile_clusters << std::endl;
        
        // Read the filter data for this tile
        uint32_t tile_passed = 0;
        for (uint32_t i = 0; i < tile_clusters; ++i) {
            char passed_char;
            filter_stream.read(&passed_char, 1);
            bool passed = (passed_char == 1);
            passing_clusters.push_back(passed);
            if (passed) {
                num_clusters_passed++;
                tile_passed++;
            }
        }
        
        std::cout << "    Tile " << filter_file.filename().string() 
                  << ": " << tile_passed << "/" << tile_clusters << " clusters passed QC" << std::endl;
        
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
    // If no clusters pass QC, use a reasonable buffer size for testing
    uint32_t buffer_size = (num_clusters_passed > 0) ? num_clusters_passed : (12 * 96); // 12 tiles * 96 clusters
    std::vector<std::vector<char>> h_bcl_data_buffers(total_cycles, std::vector<char>(buffer_size));
    std::vector<std::vector<uint8_t>> cbcl_qualities(total_cycles, std::vector<uint8_t>(buffer_size));
    
    // 4. Create a mapping of cycle number to CBCL files (handle multiple files per cycle)
    std::map<int, std::vector<fs::path>> cycle_to_cbcl_files;
    for (const auto& cbcl_path : cbcl_files) {
        // Extract cycle number from filename (e.g., C1.1/L001_1.cbcl -> cycle 1)
        std::string parent_dir = cbcl_path.parent_path().filename().string();
        if (parent_dir.starts_with("C") && parent_dir.find('.') != std::string::npos) {
            int cycle = std::stoi(parent_dir.substr(1, parent_dir.find('.') - 1));
            cycle_to_cbcl_files[cycle].push_back(cbcl_path);
        }
    }
    
    // Sort CBCL files within each cycle to ensure consistent ordering
    for (auto& [cycle, files] : cycle_to_cbcl_files) {
        std::sort(files.begin(), files.end());
    }
    
    std::cout << "Found CBCL files for " << cycle_to_cbcl_files.size() << " cycles" << std::endl;
    
    // 5. Determine actual number of clusters by processing the first cycle
    uint32_t total_clusters = 0;
    if (!cycle_to_cbcl_files.empty()) {
        // Process all CBCL files from the first cycle to get total cluster count
        const auto& first_cycle_files = cycle_to_cbcl_files.begin()->second;
        for (const auto& cbcl_path : first_cycle_files) {
            std::ifstream file(cbcl_path, std::ios::binary);
            if (file) {
                CbclHeader header = parse_cbcl_header(file);
                std::vector<CbclTileInfo> tiles = parse_cbcl_tile_list(file, header);
                for (const auto& tile : tiles) {
                    total_clusters += tile.num_clusters;
                }
            }
        }
        std::cout << "Total clusters across all tiles: " << total_clusters << std::endl;
    }
    
    // 6. Resize buffers to actual cluster count
    if (total_clusters > 0) {
        // TEMPORARY: Limit clusters for testing to avoid excessive memory usage
        const uint32_t MAX_CLUSTERS_FOR_TESTING = 10000000; // 10M clusters max
        if (total_clusters > MAX_CLUSTERS_FOR_TESTING) {
            std::cout << "WARNING: Limiting processing to first " << MAX_CLUSTERS_FOR_TESTING 
                      << " clusters out of " << total_clusters << " for testing" << std::endl;
            total_clusters = MAX_CLUSTERS_FOR_TESTING;
        }
        
        for (auto& buffer : h_bcl_data_buffers) {
            buffer.resize(total_clusters);
        }
        for (auto& buffer : cbcl_qualities) {
            buffer.resize(total_clusters);
        }
    }
    
    // 7. Process CBCL files for each cycle
    std::cout << "Reading " << cycle_to_cbcl_files.size() << " cycles using " 
              << omp_get_max_threads() << " threads for decompression..." << std::endl;
    
    for (int cycle = 1; cycle <= total_cycles; ++cycle) {
        auto it = cycle_to_cbcl_files.find(cycle);
        if (it == cycle_to_cbcl_files.end()) {
            std::cerr << "Warning: No CBCL file found for cycle " << cycle << std::endl;
            continue;
        }
        
        const auto& cbcl_paths = it->second;
        std::cout << "Processing cycle " << cycle << " with " << cbcl_paths.size() << " CBCL files" << std::endl;
        
        uint32_t cluster_offset = 0;
        
        // Process each CBCL file for this cycle
        for (const auto& cbcl_path : cbcl_paths) {
            std::cout << "  Processing file: " << cbcl_path.filename().string() << std::endl;
            
            std::ifstream cbcl_file(cbcl_path, std::ios::binary);
            if (!cbcl_file) {
                throw std::runtime_error("Failed to open CBCL file: " + cbcl_path.string());
            }
            
            try {
                // Parse CBCL header
                CbclHeader header = parse_cbcl_header(cbcl_file);
                
                // Parse tile list
                std::vector<CbclTileInfo> tiles = parse_cbcl_tile_list(cbcl_file, header);
                
                std::cout << "    Processing " << tiles.size() << " tiles from " << cbcl_path.filename().string() << std::endl;
                
                // Process each tile
                for (const auto& tile : tiles) {
                    std::cout << "      Processing tile " << tile.tile_id 
                              << " with " << tile.num_clusters << " clusters" << std::endl;
                    
                    CbclBlock block = parse_cbcl_block(cbcl_file, tile);
                    
                    if (block.basecalls.size() > 0) {
                        // Copy basecalls and qualities to the appropriate cycle buffer
                        int cycle_index = cycle - 1; // 0-based index
                        
                        for (uint32_t i = 0; i < block.basecalls.size() && cluster_offset + i < total_clusters; ++i) {
                            h_bcl_data_buffers[cycle_index][cluster_offset + i] = static_cast<char>(block.basecalls[i]);
                            
                            if (i < block.qualities.size()) {
                                cbcl_qualities[cycle_index][cluster_offset + i] = block.qualities[i];
                            }
                        }
                        
                        cluster_offset += tile.num_clusters;
                    }
                }
                
            } catch (const std::exception& e) {
                std::cerr << "Error processing CBCL file " << cbcl_path.filename().string() 
                          << ": " << e.what() << std::endl;
                throw;
            }
        }
    }

    // 8. Use CUDA to process CBCL data on GPU
    std::vector<Read> reads;
    
    // Prepare data pointers for CUDA processing
    std::vector<char*> h_bcl_data_ptrs;
    std::vector<size_t> h_bcl_sizes;
    
    for (auto& buffer : h_bcl_data_buffers) {
        h_bcl_data_ptrs.push_back(buffer.data());
        h_bcl_sizes.push_back(buffer.size());
    }
    
    std::cout << "Processing " << total_clusters << " clusters using CUDA GPU acceleration..." << std::endl;
    
    // Use CUDA to decode BCL data on GPU
    decode_bcl_data_cuda(h_bcl_data_ptrs, h_bcl_sizes, run_structure.read_segments, reads, total_clusters);
    
    std::cout << "Created " << reads.size() << " reads with sequences:" << std::endl;
    for (size_t i = 0; i < std::min(reads.size(), size_t(3)); ++i) {
        std::cout << "  Read " << i << ": R1=" << reads[i].sequence.substr(0, 10) << "..."
                  << ", I1=" << reads[i].index1 << ", I2=" << reads[i].index2 << std::endl;
    }
    
    return reads;
}
