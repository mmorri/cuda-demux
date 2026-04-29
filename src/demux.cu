#include "demux.h"
#include "bcl_parser_cuda.h"

#include <cuda_runtime.h>
#include <tinyxml2.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

#define CUDA_CHECK(expr)                                                                       \
    do {                                                                                       \
        cudaError_t _e = (expr);                                                               \
        if (_e != cudaSuccess) {                                                               \
            throw std::runtime_error(std::string("CUDA error at " __FILE__ ":") +              \
                                     std::to_string(__LINE__) + ": " + cudaGetErrorString(_e));\
        }                                                                                      \
    } while (0)

constexpr int kMaxBarcodes = 1024;
constexpr int kMaxBarcodeLen = 32;     // 2 bits/base * 32 = 64-bit packed

__constant__ uint64_t c_barcode_codes[kMaxBarcodes];
__constant__ uint64_t c_barcode_n_masks[kMaxBarcodes];
__constant__ int c_num_barcodes;
__constant__ int c_barcode_len;
__constant__ uint64_t c_barcode_mask;     // 2 bits per base, low (2*len) bits set
__constant__ uint64_t c_barcode_pair_mask;// every other bit (low 2*len-bits even positions)

bool verbose_log() {
    static const bool v = []() {
        const char* e = std::getenv("CUDA_DEMUX_VERBOSE");
        return e && e[0] && e[0] != '0';
    }();
    return v;
}

std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

}  // namespace

std::vector<SampleInfo> load_sample_info(const std::string& samplesheet) {
    std::vector<SampleInfo> samples;
    std::ifstream file(samplesheet);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open samplesheet file: " << samplesheet << std::endl;
        return samples;
    }

    std::string line, current_section;
    bool in_data_section = false;
    int sample_id_col = -1, index_col = -1, index2_col = -1, lane_col = -1;

    while (std::getline(file, line)) {
        if (!line.empty() && static_cast<unsigned char>(line[0]) == 0xEF) {
            line = line.substr(3);
        }
        line = trim(line);
        if (line.empty()) continue;

        if (line[0] == '[') {
            size_t end = line.find(']');
            if (end != std::string::npos) {
                current_section = line.substr(1, end - 1);
                in_data_section = (current_section == "BCLConvert_Data" || current_section == "Data");
                sample_id_col = index_col = index2_col = lane_col = -1;
                continue;
            }
        }

        if (!in_data_section) continue;

        std::vector<std::string> fields;
        {
            std::stringstream ss(line);
            std::string field;
            while (std::getline(ss, field, ',')) {
                fields.push_back(trim(field));
            }
        }

        if (sample_id_col == -1) {
            for (size_t i = 0; i < fields.size(); ++i) {
                std::string h = fields[i];
                std::transform(h.begin(), h.end(), h.begin(), ::tolower);
                if (h == "sample_id" || h == "sampleid" || h == "sample") sample_id_col = i;
                else if (h == "index" || h == "index1" || h == "i7_index_id") index_col = i;
                else if (h == "index2" || h == "i5_index_id") index2_col = i;
                else if (h == "lane") lane_col = i;
            }
            continue;
        }

        if (sample_id_col >= 0 && sample_id_col < static_cast<int>(fields.size())) {
            SampleInfo s;
            s.sample_id = fields[sample_id_col];
            if (index_col >= 0 && index_col < static_cast<int>(fields.size())) {
                s.index1 = fields[index_col];
            }
            if (index2_col >= 0 && index2_col < static_cast<int>(fields.size())) {
                s.index2 = fields[index2_col];
            }
            if (lane_col >= 0 && lane_col < static_cast<int>(fields.size())) {
                const std::string& ls = fields[lane_col];
                if (!ls.empty()) {
                    try { s.lane = std::stoi(ls); } catch (...) { s.lane = 0; }
                }
            }
            if (!s.sample_id.empty() && s.sample_id != "Sample_ID" && s.sample_id != "SampleID") {
                samples.push_back(std::move(s));
            }
        }
    }

    std::cout << "Loaded " << samples.size() << " samples from samplesheet" << std::endl;
    if (verbose_log()) {
        for (const auto& s : samples) {
            std::cout << "Sample: " << s.sample_id << ", Index1: " << s.index1
                      << ", Index2: " << s.index2 << std::endl;
        }
    }
    return samples;
}

bool validate_sample_barcodes(const std::vector<SampleInfo>& samples) {
    if (samples.empty()) {
        std::cerr << "Error: SampleSheet has no samples" << std::endl;
        return false;
    }
    const size_t i1 = samples.front().index1.length();
    const size_t i2 = samples.front().index2.length();
    if (i1 + i2 == 0) {
        std::cerr << "Error: SampleSheet contains empty barcode definitions" << std::endl;
        return false;
    }
    std::unordered_set<std::string> keys;
    for (const auto& s : samples) {
        if (s.index1.length() != i1 || s.index2.length() != i2) {
            std::cerr << "Error: Variable index lengths are not supported. Sample "
                      << s.sample_id << " has Index1 length " << s.index1.length()
                      << " and Index2 length " << s.index2.length()
                      << "; expected " << i1 << " and " << i2 << std::endl;
            return false;
        }
        std::string key = s.getCombinedBarcode() + "#" + std::to_string(s.lane);
        if (!keys.insert(key).second) {
            std::cerr << "Error: Duplicate barcode/lane combination in SampleSheet for "
                      << s.getCombinedBarcode() << " lane " << s.lane << std::endl;
            return false;
        }
    }
    return true;
}

namespace {

bool detect_reverse_complement_i5(const std::string& run_folder) {
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
        std::cerr << "Warning: RunParameters.xml not found at " << rp.string()
                  << ". Assuming no i5 reverse-complement." << std::endl;
        return false;
    }
    tinyxml2::XMLDocument doc;
    if (doc.LoadFile(rp.string().c_str()) != tinyxml2::XML_SUCCESS) {
        std::cerr << "Warning: Failed to parse RunParameters.xml. Assuming no i5 reverse-complement." << std::endl;
        return false;
    }
    std::string text;
    const char* keys[] = {"InstrumentName", "InstrumentType", "ApplicationName", "Platform"};
    for (const char* k : keys) {
        auto* el = doc.FirstChildElement(k);
        if (el && el->GetText()) { text = el->GetText(); break; }
    }
    if (text.empty()) {
        tinyxml2::XMLPrinter pr;
        doc.Print(&pr);
        text = pr.CStr();
    }
    auto lower = text;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    bool rc = (lower.find("nextseq") != std::string::npos) ||
              (lower.find("miniseq") != std::string::npos) ||
              (lower.find("novaseq") != std::string::npos);
    std::cout << "Instrument detection -> i5 RC = " << (rc ? "true" : "false") << std::endl;
    return rc;
}

__device__ __host__ inline uint8_t encode_base(char c, bool& is_n) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: is_n = true; return 0;
    }
}

void encode_barcode_table(const std::vector<std::string>& barcodes,
                          int barcode_len,
                          std::vector<uint64_t>& codes,
                          std::vector<uint64_t>& n_masks) {
    codes.assign(barcodes.size(), 0);
    n_masks.assign(barcodes.size(), 0);
    for (size_t b = 0; b < barcodes.size(); ++b) {
        const std::string& s = barcodes[b];
        if (static_cast<int>(s.size()) != barcode_len) {
            throw std::runtime_error("Barcode length mismatch in table");
        }
        uint64_t code = 0, nmask = 0;
        for (int i = 0; i < barcode_len; ++i) {
            bool is_n = false;
            uint8_t v = encode_base(s[i], is_n);
            code |= (static_cast<uint64_t>(v) << (2 * i));
            if (is_n) nmask |= (1ULL << (2 * i));
        }
        codes[b] = code;
        n_masks[b] = nmask;
    }
}

__global__ void pack_barcodes_kernel(const char* d_seq,
                                     int total_seq_len,
                                     int bc_offset,
                                     int bc_len,
                                     size_t batch_size,
                                     uint64_t* d_codes,
                                     uint64_t* d_nmasks) {
    size_t idx = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= batch_size) return;
    const char* p = d_seq + idx * static_cast<size_t>(total_seq_len) + bc_offset;
    uint64_t code = 0, nmask = 0;
    for (int i = 0; i < bc_len; ++i) {
        char c = p[i];
        uint8_t v = 0;
        bool is_n = false;
        switch (c) {
            case 'A': case 'a': v = 0; break;
            case 'C': case 'c': v = 1; break;
            case 'G': case 'g': v = 2; break;
            case 'T': case 't': v = 3; break;
            default: is_n = true; break;
        }
        code |= (static_cast<uint64_t>(v) << (2 * i));
        if (is_n) nmask |= (1ULL << (2 * i));
    }
    d_codes[idx] = code;
    d_nmasks[idx] = nmask;
}

__global__ void match_kernel(const uint64_t* d_codes,
                             const uint64_t* d_nmasks,
                             size_t batch_size,
                             int* d_matches) {
    size_t idx = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= batch_size) return;
    uint64_t r = d_codes[idx];
    uint64_t rn = d_nmasks[idx];
    int best = -1, best_mm = 1000, second_mm = 1000;
    for (int b = 0; b < c_num_barcodes; ++b) {
        uint64_t diff = (r ^ c_barcode_codes[b]) & c_barcode_mask;
        uint64_t lo = diff & c_barcode_pair_mask;
        uint64_t hi = (diff >> 1) & c_barcode_pair_mask;
        uint64_t pair_diff = lo | hi;
        pair_diff |= rn;                          // N in read = mismatch
        pair_diff |= c_barcode_n_masks[b];        // N in sample barcode = mismatch
        int mm = __popcll(pair_diff);
        if (mm < best_mm) {
            second_mm = best_mm;
            best_mm = mm;
            best = b;
        } else if (mm < second_mm) {
            second_mm = mm;
        }
    }
    int out = -1;
    if (best_mm <= 1 && (second_mm - best_mm) >= 1) {
        out = best;
    }
    d_matches[idx] = out;
}

uint64_t make_pair_mask(int bc_len) {
    if (bc_len <= 0) return 0;
    if (bc_len >= 32) return 0x5555555555555555ULL;
    uint64_t mask = 0;
    for (int i = 0; i < bc_len; ++i) {
        mask |= (1ULL << (2 * i));
    }
    return mask;
}

uint64_t make_full_mask(int bc_len) {
    if (bc_len <= 0) return 0;
    if (bc_len >= 32) return 0xFFFFFFFFFFFFFFFFULL;
    return (1ULL << (2 * bc_len)) - 1ULL;
}

size_t pick_batch_size(int total_seq_len, size_t lane_clusters) {
    if (lane_clusters == 0 || total_seq_len <= 0) return 0;
    size_t free_mem = 0, total_mem = 0;
    cudaError_t e = cudaMemGetInfo(&free_mem, &total_mem);
    if (e != cudaSuccess) {
        free_mem = 1ULL << 31;  // 2 GiB fallback
    }
    double frac = 0.40;
    if (const char* env = std::getenv("CUDA_DEMUX_MEM_FRACTION")) {
        try {
            double v = std::stod(env);
            if (v > 0.05 && v <= 0.95) frac = v;
        } catch (...) {}
    }
    size_t usable = static_cast<size_t>(free_mem * frac);
    // Per cluster: total_seq_len input bytes (across cycles) + 2*total_seq_len output (seq+qual)
    // + 16 bytes for packed barcodes, + 4 bytes for matches.
    size_t per_cluster = static_cast<size_t>(total_seq_len) +
                         static_cast<size_t>(total_seq_len) * 2 +
                         16 + 4;
    if (per_cluster == 0) return 0;
    size_t batch = usable / per_cluster;
    if (batch < (size_t)1 << 14) batch = (size_t)1 << 14;          // 16k floor
    if (batch > (size_t)1 << 23) batch = (size_t)1 << 23;          // 8M ceiling
    if (batch > lane_clusters) batch = lane_clusters;
    if (const char* env = std::getenv("CUDA_DEMUX_BATCH_SIZE")) {
        try {
            size_t override_v = std::stoull(env);
            if (override_v > 0) batch = std::min<size_t>(override_v, lane_clusters);
        } catch (...) {}
    }
    return batch;
}

struct PackedTable {
    std::vector<uint64_t> codes;
    std::vector<uint64_t> nmasks;
    std::vector<std::vector<int>> sample_indices_per_barcode;  // index into samples
};

PackedTable build_packed_table(const std::vector<SampleInfo>& samples,
                               bool reverse_complement_i5,
                               int& barcode_len) {
    PackedTable t;
    if (samples.empty()) {
        barcode_len = 0;
        return t;
    }
    std::vector<SampleInfo> tmp = samples;
    for (auto& s : tmp) s.reverse_complement_i2 = reverse_complement_i5;

    std::unordered_map<std::string, int> bc_to_index;
    std::vector<std::string> ordered_codes;
    for (size_t i = 0; i < tmp.size(); ++i) {
        const std::string code = tmp[i].getCombinedBarcode();
        auto it = bc_to_index.find(code);
        int slot;
        if (it == bc_to_index.end()) {
            slot = static_cast<int>(ordered_codes.size());
            ordered_codes.push_back(code);
            bc_to_index.emplace(code, slot);
            t.sample_indices_per_barcode.emplace_back();
        } else {
            slot = it->second;
        }
        t.sample_indices_per_barcode[slot].push_back(static_cast<int>(i));
    }
    barcode_len = static_cast<int>(ordered_codes.front().length());
    encode_barcode_table(ordered_codes, barcode_len, t.codes, t.nmasks);
    return t;
}

void upload_table(const PackedTable& t, int barcode_len) {
    int n = static_cast<int>(t.codes.size());
    if (n > kMaxBarcodes) {
        throw std::runtime_error("Too many distinct barcodes (max " +
                                 std::to_string(kMaxBarcodes) + ")");
    }
    if (barcode_len > kMaxBarcodeLen) {
        throw std::runtime_error("Barcode length too large (max " +
                                 std::to_string(kMaxBarcodeLen) + ")");
    }
    CUDA_CHECK(cudaMemcpyToSymbol(c_barcode_codes, t.codes.data(), n * sizeof(uint64_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_barcode_n_masks, t.nmasks.data(), n * sizeof(uint64_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_num_barcodes, &n, sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_barcode_len, &barcode_len, sizeof(int)));
    uint64_t pair_mask = make_pair_mask(barcode_len);
    uint64_t full_mask = make_full_mask(barcode_len);
    CUDA_CHECK(cudaMemcpyToSymbol(c_barcode_mask, &full_mask, sizeof(uint64_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_barcode_pair_mask, &pair_mask, sizeof(uint64_t)));
}

}  // namespace

void demux_and_write(const std::vector<LaneBclData>& lanes,
                     const std::string& samplesheet,
                     const std::string& run_folder,
                     FastqWriter& writer) {
    int device_count = 0;
    cudaError_t e = cudaGetDeviceCount(&device_count);
    if (e != cudaSuccess || device_count == 0) {
        throw std::runtime_error(std::string("No CUDA-capable GPU found: ") +
                                 cudaGetErrorString(e));
    }
    if (const char* dev_env = std::getenv("CUDA_DEMUX_DEVICE")) {
        try {
            int dev = std::stoi(dev_env);
            if (dev >= 0 && dev < device_count) cudaSetDevice(dev);
        } catch (...) {}
    }
    cudaDeviceProp prop;
    int cur_dev = 0;
    cudaGetDevice(&cur_dev);
    cudaGetDeviceProperties(&prop, cur_dev);
    std::cout << "Using GPU: " << prop.name << " with compute capability "
              << prop.major << "." << prop.minor << std::endl;

    std::vector<SampleInfo> samples = load_sample_info(samplesheet);
    if (samples.empty() || !validate_sample_barcodes(samples)) {
        throw std::runtime_error("Invalid SampleSheet");
    }

    bool rc_i5 = detect_reverse_complement_i5(run_folder);

    int bc_len_a = 0;
    PackedTable table_a = build_packed_table(samples, rc_i5, bc_len_a);
    bool try_both = false;
    if (const char* tb = std::getenv("CUDA_DEMUX_TRY_BOTH_I5")) {
        if (tb[0] && tb[0] != '0') try_both = true;
    }

    int total_seq_len = 0;
    for (const auto& l : lanes) total_seq_len = std::max(total_seq_len,
                                                          l.r1_len + l.i1_len + l.i2_len + l.r2_len);
    if (total_seq_len <= 0) {
        std::cerr << "No sequence cycles in any lane; nothing to demultiplex." << std::endl;
        return;
    }

    long long matched_total = 0;
    long long unmatched_total = 0;
    std::unordered_map<std::string, long long> per_sample_counts;

    for (const auto& lane : lanes) {
        if (lane.num_clusters == 0) continue;
        const int lane_total = lane.r1_len + lane.i1_len + lane.i2_len + lane.r2_len;
        const int bc_offset = lane.r1_len;                              // index region starts after R1
        const int bc_len = lane.i1_len + lane.i2_len;
        if (bc_len <= 0) {
            std::cerr << "Lane " << lane.lane << ": no index cycles, skipping." << std::endl;
            continue;
        }
        if (bc_len != bc_len_a) {
            throw std::runtime_error("Index length mismatch between SampleSheet and run "
                                     "(samplesheet=" + std::to_string(bc_len_a) +
                                     ", run=" + std::to_string(bc_len) + ")");
        }

        // Pick orientation: optional probe on the first batch
        bool use_flip = false;
        if (try_both) {
            int bc_len_b = 0;
            PackedTable table_b = build_packed_table(samples, !rc_i5, bc_len_b);
            // Probe the first 65k clusters with each table to decide orientation
            const size_t probe = std::min<size_t>(65536, lane.num_clusters);
            // Run a minimal probe by uploading table_a, decoding+matching probe clusters,
            // then table_b, then compare match counts.
            auto probe_score = [&](const PackedTable& t) {
                upload_table(t, bc_len_a);
                CudaDecodeContext* ctx = decode_context_create(lane);
                size_t batch = std::min(probe, pick_batch_size(lane_total, probe));
                if (batch == 0) batch = probe;
                char* d_seq = nullptr;
                char* d_qual = nullptr;
                uint64_t* d_codes = nullptr;
                uint64_t* d_nmasks = nullptr;
                int* d_matches = nullptr;
                long long score = 0;
                try {
                    CUDA_CHECK(cudaMalloc(&d_seq, batch * lane_total));
                    CUDA_CHECK(cudaMalloc(&d_qual, batch * lane_total));
                    CUDA_CHECK(cudaMalloc(&d_codes, batch * sizeof(uint64_t)));
                    CUDA_CHECK(cudaMalloc(&d_nmasks, batch * sizeof(uint64_t)));
                    CUDA_CHECK(cudaMalloc(&d_matches, batch * sizeof(int)));
                    if (!decode_bcl_batch(ctx, lane, 0, batch, d_seq, d_qual)) {
                        throw std::runtime_error("probe decode failed");
                    }
                    int threads = 256;
                    int blocks = static_cast<int>((batch + threads - 1) / threads);
                    pack_barcodes_kernel<<<blocks, threads>>>(d_seq, lane_total, bc_offset,
                                                              bc_len, batch, d_codes, d_nmasks);
                    CUDA_CHECK(cudaGetLastError());
                    match_kernel<<<blocks, threads>>>(d_codes, d_nmasks, batch, d_matches);
                    CUDA_CHECK(cudaGetLastError());
                    std::vector<int> h_matches(batch);
                    CUDA_CHECK(cudaMemcpy(h_matches.data(), d_matches, batch * sizeof(int),
                                          cudaMemcpyDeviceToHost));
                    for (int m : h_matches) if (m >= 0) ++score;
                } catch (...) {
                    if (d_seq) cudaFree(d_seq);
                    if (d_qual) cudaFree(d_qual);
                    if (d_codes) cudaFree(d_codes);
                    if (d_nmasks) cudaFree(d_nmasks);
                    if (d_matches) cudaFree(d_matches);
                    decode_context_destroy(ctx);
                    throw;
                }
                cudaFree(d_seq);
                cudaFree(d_qual);
                cudaFree(d_codes);
                cudaFree(d_nmasks);
                cudaFree(d_matches);
                decode_context_destroy(ctx);
                return score;
            };
            long long score_a = probe_score(table_a);
            long long score_b = probe_score(table_b);
            std::cout << "Lane " << lane.lane << " orientation probe: forward=" << score_a
                      << " reverse-comp=" << score_b << std::endl;
            if (score_b > score_a) use_flip = true;
        }
        PackedTable table_flip;
        if (use_flip) {
            int dummy_len = 0;
            table_flip = build_packed_table(samples, !rc_i5, dummy_len);
        }
        const PackedTable& active_table = use_flip ? table_flip : table_a;
        upload_table(active_table, bc_len_a);

        // Map barcode index -> sample id for this active table
        std::vector<std::string> bc_to_sample(active_table.codes.size(), std::string());
        for (size_t b = 0; b < active_table.sample_indices_per_barcode.size(); ++b) {
            for (int si : active_table.sample_indices_per_barcode[b]) {
                int req_lane = samples[si].lane;
                if (req_lane == 0 || req_lane == lane.lane) {
                    bc_to_sample[b] = samples[si].sample_id;
                    break;
                }
            }
        }

        size_t batch_size = pick_batch_size(lane_total, lane.num_clusters);
        std::cout << "Lane " << lane.lane << ": demuxing " << lane.num_clusters
                  << " clusters in batches of up to " << batch_size << "." << std::endl;

        // Allocate device buffers (reused across batches)
        char* d_seq = nullptr;
        char* d_qual = nullptr;
        uint64_t* d_codes = nullptr;
        uint64_t* d_nmasks = nullptr;
        int* d_matches = nullptr;
        char* h_seq = nullptr;
        char* h_qual = nullptr;
        int* h_matches = nullptr;

        auto cleanup = [&]() {
            if (d_seq) cudaFree(d_seq);
            if (d_qual) cudaFree(d_qual);
            if (d_codes) cudaFree(d_codes);
            if (d_nmasks) cudaFree(d_nmasks);
            if (d_matches) cudaFree(d_matches);
            if (h_seq) cudaFreeHost(h_seq);
            if (h_qual) cudaFreeHost(h_qual);
            if (h_matches) cudaFreeHost(h_matches);
        };

        try {
            CUDA_CHECK(cudaMalloc(&d_seq, batch_size * lane_total));
            CUDA_CHECK(cudaMalloc(&d_qual, batch_size * lane_total));
            CUDA_CHECK(cudaMalloc(&d_codes, batch_size * sizeof(uint64_t)));
            CUDA_CHECK(cudaMalloc(&d_nmasks, batch_size * sizeof(uint64_t)));
            CUDA_CHECK(cudaMalloc(&d_matches, batch_size * sizeof(int)));
            CUDA_CHECK(cudaHostAlloc((void**)&h_seq, batch_size * lane_total, cudaHostAllocDefault));
            CUDA_CHECK(cudaHostAlloc((void**)&h_qual, batch_size * lane_total, cudaHostAllocDefault));
            CUDA_CHECK(cudaHostAlloc((void**)&h_matches, batch_size * sizeof(int), cudaHostAllocDefault));
        } catch (...) {
            cleanup();
            throw;
        }

        CudaDecodeContext* ctx = decode_context_create(lane);
        try {
            const std::string undetermined = "undetermined";
            const int r1_len = lane.r1_len;
            const int r2_len = lane.r2_len;
            const int paired_offset = lane.r1_len + lane.i1_len + lane.i2_len;

            for (size_t start = 0; start < lane.num_clusters; start += batch_size) {
                const size_t this_batch = std::min(batch_size, lane.num_clusters - start);
                if (!decode_bcl_batch(ctx, lane, start, this_batch, d_seq, d_qual)) {
                    throw std::runtime_error("Decode batch allocation failed; reduce "
                                             "CUDA_DEMUX_BATCH_SIZE or MEM_FRACTION");
                }

                int threads = 256;
                int blocks = static_cast<int>((this_batch + threads - 1) / threads);
                pack_barcodes_kernel<<<blocks, threads>>>(d_seq, lane_total, bc_offset, bc_len,
                                                          this_batch, d_codes, d_nmasks);
                CUDA_CHECK(cudaGetLastError());
                match_kernel<<<blocks, threads>>>(d_codes, d_nmasks, this_batch, d_matches);
                CUDA_CHECK(cudaGetLastError());

                CUDA_CHECK(cudaMemcpy(h_seq, d_seq, this_batch * lane_total, cudaMemcpyDeviceToHost));
                CUDA_CHECK(cudaMemcpy(h_qual, d_qual, this_batch * lane_total, cudaMemcpyDeviceToHost));
                CUDA_CHECK(cudaMemcpy(h_matches, d_matches, this_batch * sizeof(int),
                                      cudaMemcpyDeviceToHost));

                for (size_t i = 0; i < this_batch; ++i) {
                    const int m = h_matches[i];
                    const std::string* sid = &undetermined;
                    if (m >= 0 && m < static_cast<int>(bc_to_sample.size()) && !bc_to_sample[m].empty()) {
                        sid = &bc_to_sample[m];
                        ++matched_total;
                        ++per_sample_counts[*sid];
                    } else {
                        ++unmatched_total;
                    }

                    const char* r1_seq = h_seq + i * lane_total;
                    const char* r1_qual = h_qual + i * lane_total;
                    const char* r2_seq = (r2_len > 0) ? (h_seq + i * lane_total + paired_offset) : nullptr;
                    const char* r2_qual = (r2_len > 0) ? (h_qual + i * lane_total + paired_offset) : nullptr;
                    writer.append(*sid, lane.lane,
                                  r1_seq, r1_qual, r1_len,
                                  r2_seq, r2_qual, r2_len);
                }

                if (verbose_log()) {
                    std::cout << "  Lane " << lane.lane << " batch [" << start << ", "
                              << (start + this_batch) << ") done." << std::endl;
                }
            }
        } catch (...) {
            decode_context_destroy(ctx);
            cleanup();
            throw;
        }
        decode_context_destroy(ctx);
        cleanup();
    }

    std::cout << "\nSample matching summary:" << std::endl;
    std::cout << "----------------------" << std::endl;
    for (const auto& s : samples) {
        auto it = per_sample_counts.find(s.sample_id);
        long long c = (it == per_sample_counts.end()) ? 0 : it->second;
        std::cout << s.sample_id << ": " << c << " reads" << std::endl;
    }
    std::cout << "Demultiplexing complete: " << matched_total << " matched, "
              << unmatched_total << " unmatched." << std::endl;
}
