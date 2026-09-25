#include "bcl_parser.h"
#include "common.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <omp.h>
#include <tinyxml2.h>
#include <zlib.h>
#ifdef CUDA_DEMUX_HAVE_LIBDEFLATE
#include <libdeflate.h>
#endif

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

void apply_run_structure(LaneBclData& lane, const RunLayout& run) {
    lane.total_cycles = run.total_cycles;
    lane.r1_len = run.r1_len;
    lane.i1_len = run.i1_len;
    lane.i2_len = run.i2_len;
    lane.r2_len = run.r2_len;
    lane.i2_reverse_complement = run.i2_reverse_complement;
    lane.read_segments = run.read_segments;
}

std::vector<fs::path> list_lane_dirs(const fs::path& basecalls_dir) {
    std::vector<fs::path> lane_dirs;
    for (const auto& entry : fs::directory_iterator(basecalls_dir)) {
        if (fs::is_directory(entry) && entry.path().filename().string().rfind("L00", 0) == 0) {
            lane_dirs.push_back(entry.path());
        }
    }
    if (lane_dirs.empty()) {
        lane_dirs.push_back(basecalls_dir);
    }
    std::sort(lane_dirs.begin(), lane_dirs.end());
    return lane_dirs;
}

int lane_number_from_dir(const fs::path& lane_dir) {
    const std::string name = lane_dir.filename().string();
    try {
        if (name.size() >= 4 && name[0] == 'L') return std::stoi(name.substr(1));
    } catch (...) {
    }
    return 1;
}

template <typename T>
void read_le(std::ifstream& file, T& out, const std::string& path, const char* what) {
    file.read(reinterpret_cast<char*>(&out), sizeof(T));
    if (file.gcount() != static_cast<std::streamsize>(sizeof(T))) {
        throw std::runtime_error("Truncated " + std::string(what) + " in " + path);
    }
}

LaneBclData parse_legacy_lane(const RunLayout& run, const fs::path& lane_dir);
LaneBclData parse_cbcl_lane(const RunLayout& run, const fs::path& lane_dir);

}  // namespace

RunLayout parse_run_layout(const std::string& bcl_folder) {
    fs::path bcl_dir(bcl_folder);
    fs::path run_info_path = bcl_dir / "RunInfo.xml";
    if (!fs::exists(run_info_path)) {
        throw std::runtime_error("RunInfo.xml not found in " + bcl_folder);
    }

    // 1. Parse RunInfo.xml to get run structure
    XMLDocument doc;
    XMLError result = doc.LoadFile(run_info_path.string().c_str());
    if (result != XML_SUCCESS) {
        throw std::runtime_error("Could not parse " + run_info_path.string() + ": " +
                                 doc.ErrorIDToName(result));
    }
    XMLElement* run_info_element = doc.FirstChildElement("RunInfo");
    XMLElement* run_element = run_info_element ? run_info_element->FirstChildElement("Run") : nullptr;
    XMLElement* reads_element = run_element ? run_element->FirstChildElement("Reads") : nullptr;
    if (!reads_element) {
        throw std::runtime_error("RunInfo.xml: <RunInfo><Run><Reads> not found");
    }

    RunLayout run;
    run.folder = bcl_folder;
    int read_count = 0;

    for (XMLElement* read_elem = reads_element->FirstChildElement("Read"); read_elem != nullptr;
         read_elem = read_elem->NextSiblingElement("Read")) {
        read_count++;

        const char* num_cycles_attr = read_elem->Attribute("NumCycles");
        const char* is_indexed_attr = read_elem->Attribute("IsIndexedRead");
        if (!num_cycles_attr || !is_indexed_attr) {
            throw std::runtime_error("RunInfo.xml: missing NumCycles/IsIndexedRead on Read " +
                                     std::to_string(read_count));
        }

        int num_cycles = std::stoi(num_cycles_attr);
        bool is_indexed = (std::string(is_indexed_attr) == "Y");
        const char* label;
        int segment_type;
        if (!is_indexed) {
            if (run.r1_len == 0) {
                run.r1_len = num_cycles; segment_type = 0; label = "Read1";
            } else if (run.r2_len == 0) {
                run.r2_len = num_cycles; segment_type = 3; label = "Read2";
            } else {
                throw std::runtime_error("RunInfo.xml: more than two non-index reads are not supported");
            }
        } else {
            if (run.i1_len == 0) {
                run.i1_len = num_cycles; segment_type = 1; label = "Index1";
            } else if (run.i2_len == 0) {
                run.i2_len = num_cycles; segment_type = 2; label = "Index2";
                // Newer RTA versions state the i5 read orientation explicitly.
                if (const char* rc = read_elem->Attribute("IsReverseComplement")) {
                    run.i2_reverse_complement = (std::string(rc) == "Y") ? 1 : 0;
                }
            } else {
                throw std::runtime_error("RunInfo.xml: more than two index reads are not supported");
            }
        }
        std::cout << "Read " << read_count << ": NumCycles=" << num_cycles
                  << ", IsIndexed=" << (is_indexed ? "Y" : "N") << " -> " << label << std::endl;

        for (int i = 0; i < num_cycles; ++i) {
            run.read_segments.push_back(segment_type);
        }
        run.total_cycles += num_cycles;
    }

    std::cout << "Run Structure: R1:" << run.r1_len << " I1:" << run.i1_len
              << " I2:" << run.i2_len << " R2:" << run.r2_len
              << " (" << run.total_cycles << " cycles)" << std::endl;
    if (run.i2_reverse_complement >= 0) {
        std::cout << "RunInfo.xml: Index2 IsReverseComplement="
                  << (run.i2_reverse_complement ? "Y" : "N") << std::endl;
    }

    // 2. Detect BCL format and list lanes
    fs::path basecalls_dir = bcl_dir / "Data" / "Intensities" / "BaseCalls";
    if (!fs::exists(basecalls_dir)) {
        throw std::runtime_error("BaseCalls directory not found: " + basecalls_dir.string());
    }

    bool has_cbcl = false;
    bool has_legacy_bcl = false;
    auto classify = [&](const fs::path& p) {
        if (p.extension() == ".cbcl") has_cbcl = true;
        else if (p.extension() == ".gz" && p.stem().extension() == ".bcl") has_legacy_bcl = true;
    };
    for (const auto& entry : fs::directory_iterator(basecalls_dir)) {
        classify(entry.path());
        const std::string lane_name = entry.path().filename().string();
        if (entry.is_directory() && lane_name.rfind("L00", 0) == 0) {
            for (const auto& sub : fs::directory_iterator(entry.path())) {
                classify(sub.path());
                const std::string cyc_name = sub.path().filename().string();
                if (sub.is_directory() && !cyc_name.empty() && cyc_name[0] == 'C') {
                    for (const auto& f : fs::directory_iterator(sub.path())) classify(f.path());
                }
            }
        }
    }

    if (has_cbcl) {
        std::cout << "CBCL format detected. Using CBCL parser." << std::endl;
        run.cbcl = true;
    } else if (has_legacy_bcl) {
        std::cout << "Legacy BCL (.bcl.gz) format detected. Using legacy BCL parser." << std::endl;
        run.cbcl = false;
    } else {
        throw std::runtime_error("No .cbcl or .bcl.gz files found under " + basecalls_dir.string());
    }
    for (const auto& d : list_lane_dirs(basecalls_dir)) run.lane_dirs.push_back(d.string());
    std::cout << "Found " << run.lane_dirs.size() << " lane(s) to process" << std::endl;
    return run;
}

LaneBclData parse_lane(const RunLayout& run, size_t lane_index) {
    if (lane_index >= run.lane_dirs.size()) {
        throw std::runtime_error("parse_lane: lane index out of range");
    }
    const fs::path lane_dir(run.lane_dirs[lane_index]);
    return run.cbcl ? parse_cbcl_lane(run, lane_dir) : parse_legacy_lane(run, lane_dir);
}

std::vector<LaneBclData> parse_bcl(const std::string& bcl_folder) {
    std::vector<LaneBclData> lanes;
    const RunLayout run = parse_run_layout(bcl_folder);
    for (size_t i = 0; i < run.lane_dirs.size(); ++i) {
        LaneBclData lane = parse_lane(run, i);
        if (lane.num_clusters > 0) lanes.push_back(std::move(lane));
    }
    return lanes;
}

// --- Legacy BCL Parser (.bcl.gz) ---
static std::vector<char> read_bcl_gz_file(const fs::path& path, uint32_t& cluster_count) {
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

namespace {

LaneBclData parse_legacy_lane(const RunLayout& run, const fs::path& lane_dir) {
    const std::string lane_name = lane_dir.filename().string();
    std::cout << "Processing lane: " << lane_name << " (legacy BCL)" << std::endl;

    LaneBclData lane;
    apply_run_structure(lane, run);
    lane.lane = lane_number_from_dir(lane_dir);

    const int total_cycles = run.total_cycles;
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
    return lane;
}

}  // namespace

// --- CBCL Parser (.cbcl) ---

CbclHeader parse_cbcl_header(std::ifstream& file, const std::string& path) {
    CbclHeader h;
    file.seekg(0);
    read_le(file, h.version, path, "CBCL version");
    read_le(file, h.header_size, path, "CBCL header size");
    read_le(file, h.bits_per_basecall, path, "CBCL bits per basecall");
    read_le(file, h.bits_per_quality, path, "CBCL bits per quality");
    if (h.bits_per_basecall != 2 || h.bits_per_quality != 2) {
        throw std::runtime_error(path + ": unsupported CBCL bit widths (basecall=" +
                                 std::to_string(h.bits_per_basecall) + ", quality=" +
                                 std::to_string(h.bits_per_quality) + "); expected 2/2");
    }

    uint32_t num_bins = 0;
    read_le(file, num_bins, path, "CBCL bin count");
    if (num_bins == 0 || num_bins > 64) {
        throw std::runtime_error(path + ": implausible CBCL quality bin count " + std::to_string(num_bins));
    }
    h.quality_bins.assign(4, 0);
    for (uint32_t i = 0; i < num_bins; ++i) {
        uint32_t from = 0, to = 0;
        read_le(file, from, path, "CBCL quality bin");
        read_le(file, to, path, "CBCL quality bin");
        if (i < 4) h.quality_bins[i] = static_cast<uint8_t>(std::min<uint32_t>(to, 63));
    }

    uint32_t num_tiles = 0;
    read_le(file, num_tiles, path, "CBCL tile count");
    h.tiles.resize(num_tiles);
    for (uint32_t i = 0; i < num_tiles; ++i) {
        CbclTileInfo& t = h.tiles[i];
        read_le(file, t.tile_id, path, "CBCL tile record");
        read_le(file, t.num_clusters, path, "CBCL tile record");
        read_le(file, t.uncompressed_block_size, path, "CBCL tile record");
        read_le(file, t.compressed_block_size, path, "CBCL tile record");
    }
    uint8_t flag = 0;
    read_le(file, flag, path, "CBCL non-PF flag");
    h.non_pf_excluded = (flag != 0);

    const auto parsed_size = static_cast<uint64_t>(file.tellg());
    if (parsed_size != h.header_size) {
        throw std::runtime_error(path + ": CBCL header size mismatch (declared " +
                                 std::to_string(h.header_size) + ", parsed " +
                                 std::to_string(parsed_size) + ")");
    }

    uint64_t offset = h.header_size;
    for (auto& t : h.tiles) {
        t.file_offset = offset;
        offset += t.compressed_block_size;
    }

    VLOG << "CBCL " << path << ": version=" << h.version << " header=" << h.header_size
         << " bins={" << int(h.quality_bins[0]) << "," << int(h.quality_bins[1]) << ","
         << int(h.quality_bins[2]) << "," << int(h.quality_bins[3]) << "} tiles=" << num_tiles
         << " non_pf_excluded=" << h.non_pf_excluded << std::endl;
    return h;
}

std::vector<uint8_t> read_cbcl_block(std::ifstream& file, const CbclTileInfo& tile,
                                     const std::string& path) {
    file.seekg(tile.file_offset);
    std::vector<char> compressed(tile.compressed_block_size);
    file.read(compressed.data(), tile.compressed_block_size);
    if (file.gcount() != static_cast<std::streamsize>(tile.compressed_block_size)) {
        throw std::runtime_error(path + ": truncated block for tile " + std::to_string(tile.tile_id));
    }

    std::vector<uint8_t> out(tile.uncompressed_block_size);
#ifdef CUDA_DEMUX_HAVE_LIBDEFLATE
    thread_local libdeflate_decompressor* dec = libdeflate_alloc_decompressor();
    if (!dec) throw std::runtime_error("libdeflate_alloc_decompressor failed");
    size_t produced = 0;
    const libdeflate_result ret = libdeflate_gzip_decompress(dec, compressed.data(), compressed.size(),
                                                             out.data(), out.size(), &produced);
    const bool ok = (ret == LIBDEFLATE_SUCCESS);
    const std::string err = "libdeflate " + std::to_string(static_cast<int>(ret));
#else
    z_stream strm = {};
    strm.avail_in = tile.compressed_block_size;
    strm.next_in = reinterpret_cast<Bytef*>(compressed.data());
    strm.avail_out = tile.uncompressed_block_size;
    strm.next_out = out.data();
    if (inflateInit2(&strm, 16 + MAX_WBITS) != Z_OK) {
        throw std::runtime_error(path + ": inflateInit2 failed");
    }
    const int ret = inflate(&strm, Z_FINISH);
    const size_t produced = strm.total_out;
    inflateEnd(&strm);
    const bool ok = (ret == Z_STREAM_END);
    const std::string err = "zlib " + std::to_string(ret);
#endif
    if (!ok || produced != tile.uncompressed_block_size) {
        throw std::runtime_error(path + ": inflate failed for tile " + std::to_string(tile.tile_id) +
                                 " (" + err + ", got " + std::to_string(produced) +
                                 " of " + std::to_string(tile.uncompressed_block_size) + " bytes)");
    }
    return out;
}

namespace {

// Where a tile's clusters live in the lane: raw index space (filter files) and
// the compact passing-filter index space used for output.
struct TileSpan {
    uint32_t raw_offset = 0;
    uint32_t raw_count = 0;
    uint32_t pf_offset = 0;
    uint32_t pf_count = 0;
};

// Parses s_<lane>_<tile>.filter -> tile id, or 0 if the name does not match.
uint32_t tile_from_filter_name(const fs::path& p, const std::string& prefix) {
    const std::string name = p.filename().string();
    if (name.rfind(prefix, 0) != 0 || p.extension() != ".filter") return 0;
    try {
        return static_cast<uint32_t>(std::stoul(name.substr(prefix.size(), name.size() - prefix.size() - 7)));
    } catch (...) {
        return 0;
    }
}

}  // namespace

namespace {

LaneBclData parse_cbcl_lane(const RunLayout& run, const fs::path& lane_dir) {
    const int total_cycles = run.total_cycles;
    {
        std::cout << "Processing lane: " << lane_dir.string() << std::endl;
        const int lane_number = lane_number_from_dir(lane_dir);

        // Per-tile filter files are named s_<lane>_<tile>.filter on every Illumina platform.
        const std::string filter_prefix = "s_" + std::to_string(lane_number) + "_";
        std::map<uint32_t, fs::path> filter_files;   // sorted by tile id
        for (const auto& entry : fs::directory_iterator(lane_dir)) {
            uint32_t tile = tile_from_filter_name(entry.path(), filter_prefix);
            if (tile) filter_files[tile] = entry.path();
        }
        if (filter_files.empty()) {
            throw std::runtime_error("No s_" + std::to_string(lane_number) +
                                     "_<tile>.filter files found in " + lane_dir.string());
        }
        std::cout << "Found " << filter_files.size() << " filter files" << std::endl;

        // pass_index_map: raw cluster index -> compact PF index (or kFiltered)
        constexpr uint32_t kFiltered = std::numeric_limits<uint32_t>::max();
        std::vector<uint32_t> pass_index_map;
        std::unordered_map<uint32_t, TileSpan> tile_spans;
        uint32_t num_clusters_total = 0;
        uint32_t num_clusters_passed = 0;

        for (const auto& [tile_id, filter_file] : filter_files) {
            std::ifstream fs_in(filter_file, std::ios::binary);
            if (!fs_in) {
                throw std::runtime_error("Failed to open filter file: " + filter_file.string());
            }
            const std::string fpath = filter_file.string();
            uint32_t first = 0, tile_clusters = 0;
            read_le(fs_in, first, fpath, "filter header");
            if (first == 0) {
                // v3 layout: uint32 0, uint32 version, uint32 clusters
                uint32_t version = 0;
                read_le(fs_in, version, fpath, "filter header");
                read_le(fs_in, tile_clusters, fpath, "filter header");
            } else {
                tile_clusters = first;
            }
            std::vector<uint8_t> flags(tile_clusters);
            fs_in.read(reinterpret_cast<char*>(flags.data()), tile_clusters);
            if (fs_in.gcount() != static_cast<std::streamsize>(tile_clusters)) {
                throw std::runtime_error("Truncated filter file: " + fpath);
            }

            TileSpan span;
            span.raw_offset = num_clusters_total;
            span.raw_count = tile_clusters;
            span.pf_offset = num_clusters_passed;
            pass_index_map.resize(num_clusters_total + tile_clusters, kFiltered);
            for (uint32_t i = 0; i < tile_clusters; ++i) {
                if (flags[i] & 1) {
                    pass_index_map[num_clusters_total + i] = num_clusters_passed++;
                    ++span.pf_count;
                }
            }
            num_clusters_total += tile_clusters;
            tile_spans[tile_id] = span;
            VLOG << "  Tile " << tile_id << ": " << span.pf_count << "/" << tile_clusters
                 << " clusters passed filter" << std::endl;
        }

        std::cout << "Filter data: " << num_clusters_total << " total clusters, "
                  << num_clusters_passed << " passing filter." << std::endl;
        LaneBclData lane;
        apply_run_structure(lane, run);
        lane.lane = lane_number;
        if (num_clusters_passed == 0) {
            std::cerr << "Warning: No clusters passed filter in lane "
                      << lane_dir.filename().string() << std::endl;
            return lane;
        }

        // cycle -> CBCL files (one per surface, sorted by name)
        std::map<int, std::vector<fs::path>> cycle_to_cbcl_files;
        for (const auto& entry : fs::directory_iterator(lane_dir)) {
            const std::string cyc_name = entry.path().filename().string();
            if (!entry.is_directory() || cyc_name.empty() || cyc_name[0] != 'C') continue;
            const size_t dot = cyc_name.find('.');
            int cycle = 0;
            try {
                cycle = std::stoi(cyc_name.substr(1, dot == std::string::npos ? std::string::npos : dot - 1));
            } catch (...) {
                continue;
            }
            for (const auto& sub : fs::directory_iterator(entry.path())) {
                if (sub.path().extension() == ".cbcl") cycle_to_cbcl_files[cycle].push_back(sub.path());
            }
        }
        for (auto& [cycle, files] : cycle_to_cbcl_files) {
            std::sort(files.begin(), files.end());
        }
        if (static_cast<int>(cycle_to_cbcl_files.size()) != total_cycles) {
            throw std::runtime_error("Lane " + lane_dir.filename().string() + ": found CBCL data for " +
                                     std::to_string(cycle_to_cbcl_files.size()) + " cycles but RunInfo.xml declares " +
                                     std::to_string(total_cycles));
        }
        for (int c = 1; c <= total_cycles; ++c) {
            if (!cycle_to_cbcl_files.count(c)) {
                throw std::runtime_error("Lane " + lane_dir.filename().string() + ": missing cycle C" +
                                         std::to_string(c) + ".1");
            }
        }

        lane.num_clusters = num_clusters_passed;
        // Allocated (and first-touched) inside the parallel loop: zero-filling
        // cycles x clusters bytes from one thread would take seconds.
        lane.bcl.resize(total_cycles);

        std::cout << "Reading " << total_cycles << " cycles using up to "
                  << omp_get_max_threads() << " threads for decompression..." << std::endl;

        std::atomic<bool> hit_error{false};
        std::string first_error;

        #pragma omp parallel for schedule(dynamic)
        for (int cycle = 1; cycle <= total_cycles; ++cycle) {
            if (hit_error.load()) continue;
            std::vector<uint8_t>& dest = lane.bcl[cycle - 1];
            try {
                dest.resize(num_clusters_passed);
                for (const auto& cbcl_path : cycle_to_cbcl_files.at(cycle)) {
                    const std::string path = cbcl_path.string();
                    std::ifstream cbcl_file(cbcl_path, std::ios::binary);
                    if (!cbcl_file) {
                        throw std::runtime_error("Failed to open CBCL file: " + path);
                    }
                    const CbclHeader header = parse_cbcl_header(cbcl_file, path);
                    const uint8_t* qmap = header.quality_bins.data();

                    for (const auto& tile : header.tiles) {
                        auto it = tile_spans.find(tile.tile_id);
                        if (it == tile_spans.end()) {
                            throw std::runtime_error(path + ": tile " + std::to_string(tile.tile_id) +
                                                     " has no filter file");
                        }
                        const TileSpan& span = it->second;
                        const uint32_t expected = header.non_pf_excluded ? span.pf_count : span.raw_count;
                        if (tile.num_clusters != expected) {
                            throw std::runtime_error(path + ": tile " + std::to_string(tile.tile_id) + " has " +
                                                     std::to_string(tile.num_clusters) + " clusters, filter file implies " +
                                                     std::to_string(expected));
                        }
                        if (tile.uncompressed_block_size < (tile.num_clusters + 1) / 2) {
                            throw std::runtime_error(path + ": tile " + std::to_string(tile.tile_id) +
                                                     " block too small for its cluster count");
                        }
                        if (tile.num_clusters == 0) continue;

                        const std::vector<uint8_t> block = read_cbcl_block(cbcl_file, tile, path);

                        // Each byte packs two clusters: low nibble first. Within a nibble,
                        // bits 0-1 are the base and bits 2-3 the quality bin. Bin 0 is a
                        // no-call and is stored as byte 0 so the decoder emits 'N'.
                        auto decode = [&](uint32_t i) -> uint8_t {
                            const uint8_t nib = (i & 1) ? (block[i >> 1] >> 4) : (block[i >> 1] & 0x0F);
                            const uint8_t q = qmap[nib >> 2];
                            return q ? static_cast<uint8_t>((q << 2) | (nib & 3)) : 0;
                        };
                        if (header.non_pf_excluded) {
                            uint8_t* out = dest.data() + span.pf_offset;
                            for (uint32_t i = 0; i < tile.num_clusters; ++i) out[i] = decode(i);
                        } else {
                            const uint32_t* pmap = pass_index_map.data() + span.raw_offset;
                            for (uint32_t i = 0; i < tile.num_clusters; ++i) {
                                const uint32_t out_idx = pmap[i];
                                if (out_idx != kFiltered) dest[out_idx] = decode(i);
                            }
                        }
                    }
                }
            } catch (const std::exception& e) {
                bool expected = false;
                if (hit_error.compare_exchange_strong(expected, true)) {
                    #pragma omp critical
                    first_error = std::string("cycle ") + std::to_string(cycle) + ": " + e.what();
                }
            }
        }

        if (hit_error.load()) {
            throw std::runtime_error("CBCL parse failed: " + first_error);
        }

        std::cout << "Lane " << lane_dir.filename().string()
                  << ": prepared " << num_clusters_passed << " clusters across "
                  << total_cycles << " cycles." << std::endl;
        return lane;
    }
}

}  // namespace
