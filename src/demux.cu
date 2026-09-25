#include "demux.h"
#include "bcl_parser_cuda.h"

#include <cuda_runtime.h>
#include <tinyxml2.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
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
constexpr int kPipelineSlots = 2;      // batches in flight: GPU(k+1) overlaps host(k)
constexpr int kMaxBarcodeLen = 32;     // 2 bits/base * 32 = 64-bit packed

__constant__ uint64_t c_barcode_codes[kMaxBarcodes];
__constant__ uint64_t c_barcode_n_masks[kMaxBarcodes];
__constant__ int c_num_barcodes;
__constant__ int c_barcode_len;
__constant__ uint64_t c_barcode_mask;     // 2 bits per base, low (2*len) bits set
__constant__ uint64_t c_barcode_pair_mask;// every other bit (low 2*len-bits even positions)
__constant__ uint64_t c_i1_pair_mask;     // pair-mask bits belonging to index 1
__constant__ uint64_t c_i2_pair_mask;     // pair-mask bits belonging to index 2
__constant__ int c_max_mm_i1;             // allowed mismatches in index 1
__constant__ int c_max_mm_i2;             // allowed mismatches in index 2

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

bool detect_reverse_complement_i5(const std::string& run_folder, int runinfo_rc) {
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
    if (runinfo_rc >= 0) {
        std::cout << "RunInfo.xml Index2 IsReverseComplement -> i5 RC = "
                  << (runinfo_rc ? "true" : "false") << std::endl;
        return runinfo_rc != 0;
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
    // Best candidate by total mismatches; ties leave the read undetermined.
    // A candidate must also respect the per-index limits (bcl-convert semantics:
    // BarcodeMismatchesIndex1 / BarcodeMismatchesIndex2 apply independently).
    int best = -1, best_mm = 1000, second_mm = 1000;
    bool best_ok = false;
    for (int b = 0; b < c_num_barcodes; ++b) {
        uint64_t diff = (r ^ c_barcode_codes[b]) & c_barcode_mask;
        uint64_t lo = diff & c_barcode_pair_mask;
        uint64_t hi = (diff >> 1) & c_barcode_pair_mask;
        uint64_t pair_diff = lo | hi;
        pair_diff |= rn;                          // N in read = mismatch
        pair_diff |= c_barcode_n_masks[b];        // N in sample barcode = mismatch
        int mm1 = __popcll(pair_diff & c_i1_pair_mask);
        int mm2 = __popcll(pair_diff & c_i2_pair_mask);
        int mm = mm1 + mm2;
        if (mm < best_mm) {
            second_mm = best_mm;
            best_mm = mm;
            best = b;
            best_ok = (mm1 <= c_max_mm_i1) && (mm2 <= c_max_mm_i2);
        } else if (mm < second_mm) {
            second_mm = mm;
        }
    }
    int out = -1;
    if (best_ok && (second_mm - best_mm) >= 1) {
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
    // Per cluster and per pipeline slot: total_seq_len input bytes (across
    // cycles) + 2*total_seq_len output (seq+qual) + 16 bytes packed barcode +
    // 4 bytes match. Two slots are live at once.
    size_t per_cluster = kPipelineSlots * (static_cast<size_t>(total_seq_len) * 3 + 16 + 4);
    if (per_cluster == 0) return 0;
    size_t batch = usable / per_cluster;
    if (batch < (size_t)1 << 14) batch = (size_t)1 << 14;          // 16k floor
    // 1M ceiling: keeps the pinned buffers (two slots of cycles + 2 x row text)
    // and the per-batch FASTQ text modest; larger batches were not faster.
    if (batch > (size_t)1 << 20) batch = (size_t)1 << 20;
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

int mismatch_limit() {
    static const int v = []() {
        int n = 1;
        if (const char* e = std::getenv("CUDA_DEMUX_MISMATCHES")) {
            try { n = std::stoi(e); } catch (...) {}
        }
        return std::max(0, std::min(n, 4));
    }();
    return v;
}

void upload_table(const PackedTable& t, int barcode_len, int i1_len) {
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
    const uint64_t i1_mask = make_pair_mask(std::min(i1_len, barcode_len));
    const uint64_t i2_mask = pair_mask & ~i1_mask;
    const int mm = mismatch_limit();
    CUDA_CHECK(cudaMemcpyToSymbol(c_i1_pair_mask, &i1_mask, sizeof(uint64_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_i2_pair_mask, &i2_mask, sizeof(uint64_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_max_mm_i1, &mm, sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_max_mm_i2, &mm, sizeof(int)));
}

}  // namespace

namespace {

// Everything one in-flight batch needs. Two of these alternate so the GPU can
// decode batch k+1 while the host formats and compresses batch k.
struct BatchBuffers {
    char* d_seq = nullptr;
    char* d_qual = nullptr;
    uint64_t* d_codes = nullptr;
    uint64_t* d_nmasks = nullptr;
    int* d_matches = nullptr;
    char* h_seq = nullptr;
    char* h_qual = nullptr;
    int* h_matches = nullptr;
    cudaEvent_t done = nullptr;
    size_t start = 0;
    size_t count = 0;

    void allocate(size_t batch, int lane_total) {
        CUDA_CHECK(cudaMalloc(&d_seq, batch * lane_total));
        CUDA_CHECK(cudaMalloc(&d_qual, batch * lane_total));
        CUDA_CHECK(cudaMalloc(&d_codes, batch * sizeof(uint64_t)));
        CUDA_CHECK(cudaMalloc(&d_nmasks, batch * sizeof(uint64_t)));
        CUDA_CHECK(cudaMalloc(&d_matches, batch * sizeof(int)));
        CUDA_CHECK(cudaHostAlloc((void**)&h_seq, batch * lane_total, cudaHostAllocDefault));
        CUDA_CHECK(cudaHostAlloc((void**)&h_qual, batch * lane_total, cudaHostAllocDefault));
        CUDA_CHECK(cudaHostAlloc((void**)&h_matches, batch * sizeof(int), cudaHostAllocDefault));
        CUDA_CHECK(cudaEventCreateWithFlags(&done, cudaEventDisableTiming));
    }
    void release() {
        if (d_seq) cudaFree(d_seq);
        if (d_qual) cudaFree(d_qual);
        if (d_codes) cudaFree(d_codes);
        if (d_nmasks) cudaFree(d_nmasks);
        if (d_matches) cudaFree(d_matches);
        if (h_seq) cudaFreeHost(h_seq);
        if (h_qual) cudaFreeHost(h_qual);
        if (h_matches) cudaFreeHost(h_matches);
        if (done) cudaEventDestroy(done);
        *this = BatchBuffers();
    }
};

// Enqueues decode + barcode pack + match + copy-back for one batch on the
// context stream and records the slot's completion event.
void enqueue_batch(CudaDecodeContext* ctx, const LaneBclData& lane, BatchBuffers& b,
                   size_t start, size_t count, int slot, int lane_total, int bc_offset, int bc_len) {
    b.start = start;
    b.count = count;
    if (!decode_bcl_batch(ctx, lane, start, count, slot, b.d_seq, b.d_qual)) {
        throw std::runtime_error("Decode batch allocation failed; reduce "
                                 "CUDA_DEMUX_BATCH_SIZE or MEM_FRACTION");
    }
    cudaStream_t stream = decode_context_stream(ctx);
    const int threads = 256;
    const int blocks = static_cast<int>((count + threads - 1) / threads);
    pack_barcodes_kernel<<<blocks, threads, 0, stream>>>(b.d_seq, lane_total, bc_offset, bc_len,
                                                         count, b.d_codes, b.d_nmasks);
    CUDA_CHECK(cudaGetLastError());
    match_kernel<<<blocks, threads, 0, stream>>>(b.d_codes, b.d_nmasks, count, b.d_matches);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(b.h_seq, b.d_seq, count * lane_total, cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(b.h_qual, b.d_qual, count * lane_total, cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(b.h_matches, b.d_matches, count * sizeof(int), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaEventRecord(b.done, stream));
}

// Counts matches of `table` on the first `probe` clusters of the lane.
long long probe_orientation_score(const PackedTable& table, const LaneBclData& lane,
                                  size_t probe, int lane_total, int bc_offset, int bc_len) {
    upload_table(table, bc_len, lane.i1_len);
    BatchBuffers b;
    CudaDecodeContext* ctx = nullptr;
    long long score = 0;
    try {
        b.allocate(probe, lane_total);
        ctx = decode_context_create(lane.total_cycles, probe, 1);
        decode_context_set_lane(ctx, lane);
        enqueue_batch(ctx, lane, b, 0, probe, 0, lane_total, bc_offset, bc_len);
        CUDA_CHECK(cudaEventSynchronize(b.done));
        for (size_t i = 0; i < probe; ++i) if (b.h_matches[i] >= 0) ++score;
    } catch (...) {
        decode_context_destroy(ctx);
        b.release();
        throw;
    }
    decode_context_destroy(ctx);
    b.release();
    return score;
}

}  // namespace

struct Demuxer::Impl {
    FastqWriter& writer;
    std::vector<SampleInfo> samples;
    bool rc_i5 = false;
    bool try_both = false;
    int bc_len_a = 0;
    PackedTable table_a;
    int total_seq_len = 0;
    int num_cycles = 0;
    size_t capacity = 0;
    BatchBuffers bufs[kPipelineSlots];
    CudaDecodeContext* ctx = nullptr;
    long long matched_total = 0;
    long long unmatched_total = 0;
    std::unordered_map<std::string, long long> per_sample_counts;

    Impl(const RunLayout& run, const std::string& samplesheet, FastqWriter& w) : writer(w) {
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

        samples = load_sample_info(samplesheet);
        if (samples.empty() || !validate_sample_barcodes(samples)) {
            throw std::runtime_error("Invalid SampleSheet");
        }

        rc_i5 = detect_reverse_complement_i5(run.folder, run.i2_reverse_complement);
        std::cout << "Barcode mismatches allowed per index: " << mismatch_limit() << std::endl;
        table_a = build_packed_table(samples, rc_i5, bc_len_a);
        if (const char* tb = std::getenv("CUDA_DEMUX_TRY_BOTH_I5")) {
            if (tb[0] && tb[0] != '0') try_both = true;
        }

        total_seq_len = run.r1_len + run.i1_len + run.i2_len + run.r2_len;
        num_cycles = run.total_cycles;
        if (total_seq_len <= 0) {
            throw std::runtime_error("No sequence cycles in run; nothing to demultiplex.");
        }
        const int bc_len = run.i1_len + run.i2_len;
        if (bc_len <= 0) {
            throw std::runtime_error("Run has no index cycles; nothing to demultiplex.");
        }
        if (bc_len != bc_len_a) {
            throw std::runtime_error("Index length mismatch between SampleSheet and run "
                                     "(samplesheet=" + std::to_string(bc_len_a) +
                                     ", run=" + std::to_string(bc_len) + ")");
        }

        // Pipeline buffers are sized once and reused by every lane: pinning
        // several GB of host memory costs seconds per allocation.
        capacity = pick_batch_size(total_seq_len, std::numeric_limits<size_t>::max());
        try {
            for (auto& b : bufs) b.allocate(capacity, total_seq_len);
            ctx = decode_context_create(num_cycles, capacity, kPipelineSlots);
        } catch (...) {
            release();
            throw;
        }
    }

    ~Impl() { release(); }

    void release() {
        decode_context_destroy(ctx);
        ctx = nullptr;
        for (auto& b : bufs) b.release();
    }

    void process_lane(const LaneBclData& lane) {
        if (lane.num_clusters == 0) return;
        const int lane_total = lane.r1_len + lane.i1_len + lane.i2_len + lane.r2_len;
        const int bc_offset = lane.r1_len;                              // index region starts after R1
        const int bc_len = lane.i1_len + lane.i2_len;
        if (lane_total != total_seq_len || lane.total_cycles != num_cycles) {
            throw std::runtime_error("Lane " + std::to_string(lane.lane) +
                                     " does not match the run structure");
        }

        // Pick orientation: optional probe on the first clusters with both tables.
        bool use_flip = false;
        if (try_both) {
            int bc_len_b = 0;
            PackedTable table_b = build_packed_table(samples, !rc_i5, bc_len_b);
            const size_t probe = std::min<size_t>(65536, lane.num_clusters);
            long long score_a = probe_orientation_score(table_a, lane, probe, lane_total, bc_offset, bc_len);
            long long score_b = probe_orientation_score(table_b, lane, probe, lane_total, bc_offset, bc_len);
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
        upload_table(active_table, bc_len_a, lane.i1_len);

        // Output slots for this lane: one per distinct sample id (samplesheet
        // order, restricted to samples assigned to this lane) plus "undetermined".
        std::vector<std::string> slot_ids;
        std::unordered_map<std::string, int> slot_of_sample;
        for (const auto& s : samples) {
            if (s.lane != 0 && s.lane != lane.lane) continue;
            if (slot_of_sample.emplace(s.sample_id, static_cast<int>(slot_ids.size())).second) {
                slot_ids.push_back(s.sample_id);
            }
        }
        const int undetermined_slot = static_cast<int>(slot_ids.size());
        slot_ids.push_back("undetermined");
        const bool paired = lane.r2_len > 0;
        for (const auto& id : slot_ids) writer.ensure_open(id, lane.lane, paired);

        // barcode index -> output slot
        std::vector<int> bc_to_slot(active_table.codes.size(), undetermined_slot);
        for (size_t b = 0; b < active_table.sample_indices_per_barcode.size(); ++b) {
            for (int si : active_table.sample_indices_per_barcode[b]) {
                int req_lane = samples[si].lane;
                if (req_lane == 0 || req_lane == lane.lane) {
                    bc_to_slot[b] = slot_of_sample.at(samples[si].sample_id);
                    break;
                }
            }
        }
        std::vector<long long> slot_counts(slot_ids.size(), 0);

        const size_t batch_size = std::min(capacity, lane.num_clusters);
        const size_t num_batches = (lane.num_clusters + batch_size - 1) / batch_size;
        std::cout << "Lane " << lane.lane << ": demuxing " << lane.num_clusters
                  << " clusters in " << num_batches << " batches of up to " << batch_size << "." << std::endl;

        std::vector<int> h_slots(batch_size);
        double t_wait = 0.0, t_write = 0.0, t_enqueue = 0.0;
        using clock = std::chrono::steady_clock;

        decode_context_set_lane(ctx, lane);
        const int paired_offset = lane.r1_len + lane.i1_len + lane.i2_len;
        auto enqueue = [&](size_t k) {
            const size_t start = k * batch_size;
            const size_t count = std::min(batch_size, lane.num_clusters - start);
            enqueue_batch(ctx, lane, bufs[k % kPipelineSlots], start, count,
                          static_cast<int>(k % kPipelineSlots), lane_total, bc_offset, bc_len);
        };

        enqueue(0);
        for (size_t k = 0; k < num_batches; ++k) {
            // Keep the GPU busy with the next batch while this one is written.
            const auto t_e = clock::now();
            if (k + 1 < num_batches) enqueue(k + 1);
            t_enqueue += std::chrono::duration<double>(clock::now() - t_e).count();

            BatchBuffers& b = bufs[k % kPipelineSlots];
            const auto t0 = clock::now();
            CUDA_CHECK(cudaEventSynchronize(b.done));
            const auto t1 = clock::now();

            for (size_t i = 0; i < b.count; ++i) {
                const int m = b.h_matches[i];
                const int slot = (m >= 0 && m < static_cast<int>(bc_to_slot.size()))
                                     ? bc_to_slot[m] : undetermined_slot;
                h_slots[i] = slot;
                ++slot_counts[slot];
            }

            FastqBatch fb;
            fb.lane = lane.lane;
            fb.count = b.count;
            fb.seq = b.h_seq;
            fb.qual = b.h_qual;
            fb.stride = static_cast<size_t>(lane_total);
            fb.r1_offset = 0;
            fb.r1_len = lane.r1_len;
            fb.r2_offset = paired_offset;
            fb.r2_len = lane.r2_len;
            fb.sample_slot = h_slots.data();
            fb.sample_ids = &slot_ids;
            writer.append_batch(fb);

            t_wait += std::chrono::duration<double>(t1 - t0).count();
            t_write += std::chrono::duration<double>(clock::now() - t1).count();
            if (verbose_log()) {
                std::cout << "  Lane " << lane.lane << " batch [" << b.start << ", "
                          << (b.start + b.count) << ") done." << std::endl;
            }
        }
        if (verbose_log()) {
            std::cout << std::fixed << std::setprecision(2) << "  Lane " << lane.lane
                      << " timing: staging+enqueue " << t_enqueue << " s, GPU wait " << t_wait
                      << " s, FASTQ format+compress+write " << t_write << " s" << std::endl;
        }

        for (size_t sl = 0; sl < slot_ids.size(); ++sl) {
            if (static_cast<int>(sl) == undetermined_slot) {
                unmatched_total += slot_counts[sl];
            } else {
                matched_total += slot_counts[sl];
                per_sample_counts[slot_ids[sl]] += slot_counts[sl];
            }
        }
        std::cout << "Lane " << lane.lane << ": " << (lane.num_clusters - slot_counts[undetermined_slot])
                  << " of " << lane.num_clusters << " reads assigned ("
                  << std::fixed << std::setprecision(2)
                  << (100.0 * (lane.num_clusters - slot_counts[undetermined_slot]) / lane.num_clusters)
                  << "%)" << std::endl;
    }

    void print_summary() const {
        std::cout << "\nSample matching summary:" << std::endl;
        std::cout << "----------------------" << std::endl;
        std::unordered_set<std::string> printed;
        for (const auto& s : samples) {
            if (!printed.insert(s.sample_id).second) continue;
            auto it = per_sample_counts.find(s.sample_id);
            long long c = (it == per_sample_counts.end()) ? 0 : it->second;
            std::cout << s.sample_id << ": " << c << " reads" << std::endl;
        }
        std::cout << "Demultiplexing complete: " << matched_total << " matched, "
                  << unmatched_total << " unmatched." << std::endl;
    }
};

Demuxer::Demuxer(const RunLayout& run, const std::string& samplesheet, FastqWriter& writer)
    : impl_(std::make_unique<Impl>(run, samplesheet, writer)) {}
Demuxer::~Demuxer() = default;
void Demuxer::process_lane(const LaneBclData& lane) { impl_->process_lane(lane); }
void Demuxer::print_summary() const { impl_->print_summary(); }
