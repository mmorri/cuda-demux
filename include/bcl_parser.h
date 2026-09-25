#ifndef BCL_PARSER_H
#define BCL_PARSER_H

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include "common.h"

// One tile record from a CBCL header.
struct CbclTileInfo {
    uint32_t tile_id;
    uint32_t num_clusters;             // raw clusters, or PF clusters when non_pf_excluded
    uint32_t uncompressed_block_size;
    uint32_t compressed_block_size;
    uint64_t file_offset;              // absolute offset of the gzip block
};

// CBCL header (RTA3 format):
//   uint16 version, uint32 header_size, uint8 bits_per_basecall, uint8 bits_per_quality,
//   uint32 num_bins, num_bins x {uint32 from, uint32 to}, uint32 num_tiles,
//   num_tiles x {uint32 tile, uint32 clusters, uint32 uncompressed, uint32 compressed},
//   uint8 non_pf_clusters_excluded
struct CbclHeader {
    uint16_t version = 0;
    uint32_t header_size = 0;
    uint8_t bits_per_basecall = 0;
    uint8_t bits_per_quality = 0;
    std::vector<uint8_t> quality_bins;   // bin index -> Q-score
    std::vector<CbclTileInfo> tiles;
    bool non_pf_excluded = false;
};

// Parses the whole CBCL header (including the tile list) from an open stream.
CbclHeader parse_cbcl_header(std::ifstream& file, const std::string& path);

// Reads and inflates one tile block. Returns the raw nibble-packed bytes
// (two clusters per byte: low nibble first; bits 0-1 base, bits 2-3 q-bin).
std::vector<uint8_t> read_cbcl_block(std::ifstream& file, const CbclTileInfo& tile,
                                     const std::string& path);

// Run-level layout from RunInfo.xml plus the discovered lane directories.
struct RunLayout {
    std::string folder;
    int total_cycles = 0;
    int r1_len = 0;
    int i1_len = 0;
    int i2_len = 0;
    int r2_len = 0;
    int i2_reverse_complement = -1;   // RunInfo IsReverseComplement on Index2: -1 unknown
    std::vector<int> read_segments;   // per cycle: 0=R1, 1=I1, 2=I2, 3=R2
    bool cbcl = true;                 // false -> legacy .bcl.gz
    std::vector<std::string> lane_dirs;
};

// Reads RunInfo.xml, detects the BCL flavour and lists lane directories. Throws.
RunLayout parse_run_layout(const std::string& folder);

// Loads one lane (all cycles, passing-filter clusters only) into memory. Throws.
// Lanes are loaded one at a time so peak memory is one lane, not the run.
LaneBclData parse_lane(const RunLayout& run, size_t lane_index);

// Convenience: every lane of the run at once.
std::vector<LaneBclData> parse_bcl(const std::string& folder);

#endif
