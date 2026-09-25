#include "fastq_writer.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <utility>

#include <omp.h>
#include <zlib.h>
#ifdef CUDA_DEMUX_HAVE_LIBDEFLATE
#include <libdeflate.h>
#endif

namespace fs = std::filesystem;

namespace {
constexpr size_t kWriteBuf = 1 << 20;
// Reads per formatting/compression task. ~32k reads of 2x150 is ~11 MB raw,
// large enough for deflate to be efficient and small enough to spread across cores.
constexpr size_t kChunkReads = 32768;

std::string lane_tag(int lane) {
    char buf[16];
    std::snprintf(buf, sizeof(buf), "L%03d", lane);
    return std::string(buf);
}

std::string stream_key(const std::string& sample_id, int lane) {
    // Explicit separator: "Sample1"/lane 1 must not collide with "Sample"/lane 11.
    return sample_id + '\x1f' + std::to_string(lane);
}

void append_number(std::string& out, long long n) {
    char buf[24];
    int len = std::snprintf(buf, sizeof(buf), "%lld", n);
    out.append(buf, len);
}

// Formats reads [begin, end) of `order` for one mate into FASTQ text.
void format_records(const FastqBatch& b, const std::vector<uint32_t>& order,
                    size_t begin, size_t end, long long first_ordinal,
                    const std::string& sample_id, int mate, std::string& out) {
    const int offset = (mate == 2) ? b.r2_offset : b.r1_offset;
    const int len = (mate == 2) ? b.r2_len : b.r1_len;
    const bool paired = b.r2_len > 0;
    const size_t per_read = 1 + sample_id.size() + 1 + 20 + 3 + len + 3 + len + 1;
    out.clear();
    out.reserve((end - begin) * per_read);
    long long n = first_ordinal;
    for (size_t k = begin; k < end; ++k, ++n) {
        const size_t row = static_cast<size_t>(order[k]) * b.stride;
        out.push_back('@');
        out.append(sample_id);
        out.push_back('_');
        append_number(out, n);
        if (paired) {
            out.push_back('/');
            out.push_back(mate == 2 ? '2' : '1');
        }
        out.push_back('\n');
        out.append(b.seq + row + offset, len);
        out.append("\n+\n", 3);
        out.append(b.qual + row + offset, len);
        out.push_back('\n');
    }
}

// Compresses `in` into a standalone gzip member.
#ifdef CUDA_DEMUX_HAVE_LIBDEFLATE
void gzip_member(const std::string& in, int level, std::string& out) {
    // One compressor per thread; libdeflate compressors are not thread-safe
    // but are cheap to keep around.
    thread_local libdeflate_compressor* comp = nullptr;
    thread_local int comp_level = -1;
    if (!comp || comp_level != level) {
        if (comp) libdeflate_free_compressor(comp);
        comp = libdeflate_alloc_compressor(level);
        comp_level = level;
        if (!comp) throw std::runtime_error("libdeflate_alloc_compressor failed");
    }
    const size_t bound = libdeflate_gzip_compress_bound(comp, in.size());
    if (out.size() < bound) out.resize(bound);
    const size_t produced = libdeflate_gzip_compress(comp, in.data(), in.size(), &out[0], bound);
    if (produced == 0) throw std::runtime_error("libdeflate_gzip_compress failed");
    out.resize(produced);
}
#else
void gzip_member(const std::string& in, int level, std::string& out) {
    z_stream strm = {};
    if (deflateInit2(&strm, level, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        throw std::runtime_error("deflateInit2 failed");
    }
    const uLong bound = deflateBound(&strm, static_cast<uLong>(in.size()));
    if (out.size() < bound) out.resize(bound);
    strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(in.data()));
    strm.avail_in = static_cast<uInt>(in.size());
    strm.next_out = reinterpret_cast<Bytef*>(&out[0]);
    strm.avail_out = static_cast<uInt>(bound);
    int rc = deflate(&strm, Z_FINISH);
    const uLong produced = strm.total_out;
    deflateEnd(&strm);
    if (rc != Z_STREAM_END) {
        throw std::runtime_error("deflate failed (zlib " + std::to_string(rc) + ")");
    }
    out.resize(produced);
}
#endif

struct Task {
    int slot;
    int mate;
    size_t begin;
    size_t end;
    long long first_ordinal;
};
}  // namespace

FastqWriter::FastqWriter(std::string output_folder, bool gzip_output, int gzip_level)
    : folder_(std::move(output_folder)), gzip_(gzip_output),
      gzip_level_(std::clamp(gzip_level, 1, 9)) {
    fs::create_directories(folder_);
}

FastqWriter::~FastqWriter() {
    try {
        close();
    } catch (...) {
        // best-effort during destruction; errors during close are reported by close()
    }
}

FastqWriter::LaneStream& FastqWriter::get_or_open(const std::string& sample_id,
                                                  int lane,
                                                  bool paired) {
    const std::string key = stream_key(sample_id, lane);
    auto it = streams_.find(key);
    if (it != streams_.end()) {
        return it->second;
    }

    LaneStream s;
    s.paired = paired;
    const std::string base = folder_ + "/" + sample_id + "_" + lane_tag(lane);
    const std::string ext = gzip_ ? ".fastq.gz" : ".fastq";

    auto open_one = [&](const std::string& path) {
        std::FILE* fp = std::fopen(path.c_str(), "wb");
        if (!fp) {
            throw std::runtime_error("Could not open output file: " + path +
                                     " (" + std::strerror(errno) + ")");
        }
        std::setvbuf(fp, nullptr, _IOFBF, kWriteBuf);
        return fp;
    };

    s.r1 = open_one(base + "_R1_001" + ext);
    if (paired) {
        s.r2 = open_one(base + "_R2_001" + ext);
    }

    distinct_samples_.insert(sample_id);
    auto [ins, _] = streams_.emplace(key, std::move(s));
    return ins->second;
}

void FastqWriter::ensure_open(const std::string& sample_id, int lane, bool paired) {
    get_or_open(sample_id, lane, paired);
}

void FastqWriter::write_all(std::FILE* fp, const std::string& data, const std::string& what) {
    if (data.empty()) return;
    if (std::fwrite(data.data(), 1, data.size(), fp) != data.size()) {
        throw std::runtime_error("write failed for " + what + ": " + std::strerror(errno));
    }
}

void FastqWriter::append_batch(const FastqBatch& b) {
    if (b.count == 0) return;
    if (!b.sample_ids || !b.sample_slot || !b.seq || !b.qual) {
        throw std::runtime_error("append_batch: incomplete batch description");
    }
    const int num_slots = static_cast<int>(b.sample_ids->size());
    const bool paired = b.r2_len > 0;

    // Counting sort of reads by sample slot, preserving cluster order.
    std::vector<size_t> counts(num_slots + 1, 0);
    for (size_t i = 0; i < b.count; ++i) {
        const int s = b.sample_slot[i];
        if (s < 0 || s >= num_slots) {
            throw std::runtime_error("append_batch: sample slot out of range");
        }
        ++counts[s + 1];
    }
    for (int s = 0; s < num_slots; ++s) counts[s + 1] += counts[s];
    std::vector<uint32_t> order(b.count);
    {
        std::vector<size_t> next(counts.begin(), counts.end() - 1);
        for (size_t i = 0; i < b.count; ++i) {
            order[next[b.sample_slot[i]]++] = static_cast<uint32_t>(i);
        }
    }

    // Open streams (serial: mutates the map) and build the task list.
    std::vector<LaneStream*> streams(num_slots, nullptr);
    std::vector<Task> tasks;
    for (int s = 0; s < num_slots; ++s) {
        const size_t begin = counts[s], end = counts[s + 1];
        if (begin == end) continue;
        streams[s] = &get_or_open((*b.sample_ids)[s], b.lane, paired);
        for (size_t c = begin; c < end; c += kChunkReads) {
            const size_t ce = std::min(end, c + kChunkReads);
            const long long ordinal = streams[s]->read_count + static_cast<long long>(c - begin);
            tasks.push_back(Task{s, 1, c, ce, ordinal});
            if (paired) tasks.push_back(Task{s, 2, c, ce, ordinal});
        }
    }

    // Format (and compress) every chunk in parallel.
    if (chunk_pool_.size() < tasks.size()) chunk_pool_.resize(tasks.size());
    if (format_pool_.size() < static_cast<size_t>(omp_get_max_threads())) {
        format_pool_.resize(omp_get_max_threads());
    }
    std::string error;
    #pragma omp parallel
    {
        std::string& raw = format_pool_[omp_get_thread_num()];
        #pragma omp for schedule(dynamic)
        for (int t = 0; t < static_cast<int>(tasks.size()); ++t) {
            const Task& task = tasks[t];
            std::string& out = chunk_pool_[t];
            try {
                if (gzip_) {
                    format_records(b, order, task.begin, task.end, task.first_ordinal,
                                   (*b.sample_ids)[task.slot], task.mate, raw);
                    gzip_member(raw, gzip_level_, out);
                } else {
                    format_records(b, order, task.begin, task.end, task.first_ordinal,
                                   (*b.sample_ids)[task.slot], task.mate, out);
                }
            } catch (const std::exception& e) {
                #pragma omp critical
                if (error.empty()) error = e.what();
            }
        }
    }
    if (!error.empty()) {
        throw std::runtime_error("append_batch: " + error);
    }

    // Sequential I/O in slot/cluster order.
    for (size_t t = 0; t < tasks.size(); ++t) {
        LaneStream& s = *streams[tasks[t].slot];
        write_all(tasks[t].mate == 2 ? s.r2 : s.r1, chunk_pool_[t], (*b.sample_ids)[tasks[t].slot]);
    }
    for (int s = 0; s < num_slots; ++s) {
        if (streams[s]) streams[s]->read_count += static_cast<long long>(counts[s + 1] - counts[s]);
    }
}

void FastqWriter::append(const std::string& sample_id, int lane,
                         const char* r1_seq, const char* r1_qual, int r1_len,
                         const char* r2_seq, const char* r2_qual, int r2_len) {
    const bool paired = (r2_len > 0 && r2_seq != nullptr);
    if (!paired) r2_len = 0;
    std::string seq(r1_seq, r1_len), qual(r1_qual, r1_len);
    if (paired) {
        seq.append(r2_seq, r2_len);
        qual.append(r2_qual, r2_len);
    }
    const std::vector<std::string> ids{sample_id};
    const int slot = 0;
    FastqBatch b;
    b.lane = lane;
    b.count = 1;
    b.seq = seq.data();
    b.qual = qual.data();
    b.stride = seq.size();
    b.r1_offset = 0;
    b.r1_len = r1_len;
    b.r2_offset = r1_len;
    b.r2_len = r2_len;
    b.sample_slot = &slot;
    b.sample_ids = &ids;
    append_batch(b);
}

void FastqWriter::close() {
    if (closed_) {
        return;
    }
    closed_ = true;

    std::string error_msg;
    std::string empty_member;
    if (gzip_) gzip_member(std::string(), gzip_level_, empty_member);

    auto close_one = [&](std::FILE*& fp, long long reads) {
        if (!fp) return;
        // A zero-read gzip output must still be a valid (empty) gzip stream.
        if (gzip_ && reads == 0 && std::fwrite(empty_member.data(), 1, empty_member.size(), fp)
                                       != empty_member.size() && error_msg.empty()) {
            error_msg = std::string("write failed: ") + std::strerror(errno);
        }
        if (std::fflush(fp) != 0 && error_msg.empty()) {
            error_msg = std::string("fflush failed: ") + std::strerror(errno);
        }
        if (std::fclose(fp) != 0 && error_msg.empty()) {
            error_msg = std::string("fclose failed: ") + std::strerror(errno);
        }
        fp = nullptr;
    };
    for (auto& [_, s] : streams_) {
        close_one(s.r1, s.read_count);
        close_one(s.r2, s.read_count);
    }
    streams_.clear();

    if (!error_msg.empty()) {
        throw std::runtime_error(error_msg);
    }
}
