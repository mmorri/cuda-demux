#ifndef BCL_PARSER_H
#define BCL_PARSER_H

#include <vector>
#include <cstdint>
#include "common.h"

// CBCL-specific structures based on Picard's implementation
struct CbclHeader {
    uint32_t version;
    uint32_t header_size;
    uint32_t compression_type;
    uint32_t num_cycles;
    uint32_t num_tiles;
    uint32_t num_clusters_per_tile;
    uint32_t bits_per_basecall;
    uint32_t bits_per_quality;
    uint32_t bits_per_filter;
    uint32_t num_bases_per_cycle;
    uint32_t num_quality_bins;
    uint32_t num_filter_bins;
    uint32_t num_tiles_per_cycle;
    uint32_t tile_list_offset;
    uint32_t tile_list_size;
    uint32_t data_offset;
    uint32_t data_size;
};

struct CbclTileInfo {
    uint32_t tile_id;
    uint32_t num_clusters;
    uint32_t uncompressed_block_size;
    uint32_t compressed_block_size;
    uint64_t file_offset;
};

struct CbclBlock {
    uint32_t num_clusters;
    uint32_t num_bases;
    uint32_t num_qualities;
    uint32_t num_filters;
    std::vector<uint8_t> basecalls;
    std::vector<uint8_t> qualities;
    std::vector<uint8_t> filters;
};

std::vector<Read> parse_bcl(const std::string& folder);

// CBCL parsing functions
CbclHeader parse_cbcl_header(std::ifstream& file);
std::vector<CbclTileInfo> parse_cbcl_tile_list(std::ifstream& file, const CbclHeader& header);
CbclBlock parse_cbcl_block(std::ifstream& file, const CbclTileInfo& tile_info);
std::vector<uint8_t> decompress_cbcl_data(const std::vector<char>& compressed_data, uint32_t uncompressed_size);

#endif // BCL_PARSER_H