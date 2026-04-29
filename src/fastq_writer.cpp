#include "fastq_writer.h"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <utility>

namespace fs = std::filesystem;

namespace {
constexpr int kGzCompressionLevel = 1;
constexpr size_t kGzWriteBuf = 1 << 20;
constexpr size_t kPlainWriteBuf = 1 << 20;

std::string lane_tag(int lane) {
    char buf[16];
    std::snprintf(buf, sizeof(buf), "L%03d", lane);
    return std::string(buf);
}

std::string stream_key(const std::string& sample_id, int lane) {
    return sample_id + "\0" + std::to_string(lane);
}
}

FastqWriter::FastqWriter(std::string output_folder, bool gzip_output)
    : folder_(std::move(output_folder)), gzip_(gzip_output) {
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

    auto open_one = [&](const std::string& path, gzFile& gz, std::FILE*& fp) {
        if (gzip_) {
            gz = gzopen(path.c_str(), "wb1");
            if (!gz) {
                throw std::runtime_error("Could not open gzip output file: " + path +
                                         " (" + std::strerror(errno) + ")");
            }
            gzbuffer(gz, kGzWriteBuf);
        } else {
            fp = std::fopen(path.c_str(), "wb");
            if (!fp) {
                throw std::runtime_error("Could not open output file: " + path +
                                         " (" + std::strerror(errno) + ")");
            }
            std::setvbuf(fp, nullptr, _IOFBF, kPlainWriteBuf);
        }
    };

    open_one(base + "_R1_001" + ext, s.r1_gz, s.r1_plain);
    if (paired) {
        open_one(base + "_R2_001" + ext, s.r2_gz, s.r2_plain);
    }

    distinct_samples_.insert(sample_id);
    auto [ins, _] = streams_.emplace(key, std::move(s));
    return ins->second;
}

void FastqWriter::write_record_gz(gzFile gz, const std::string& sample_id, long long n,
                                  char suffix, const char* seq, int seq_len,
                                  const char* qual, int qual_len) {
    char header[256];
    int header_len = 0;
    if (suffix) {
        header_len = std::snprintf(header, sizeof(header), "@%s_%lld/%c\n",
                                   sample_id.c_str(), n, suffix);
    } else {
        header_len = std::snprintf(header, sizeof(header), "@%s_%lld\n",
                                   sample_id.c_str(), n);
    }
    if (header_len < 0 || header_len >= static_cast<int>(sizeof(header))) {
        throw std::runtime_error("FASTQ header overflow for sample " + sample_id);
    }

    auto must = [&](int written, int expected, const char* label) {
        if (written != expected) {
            throw std::runtime_error(std::string("gzwrite failed for ") + label +
                                     " on sample " + sample_id);
        }
    };

    must(gzwrite(gz, header, header_len), header_len, "header");
    must(gzwrite(gz, seq, seq_len), seq_len, "sequence");
    must(gzwrite(gz, "\n+\n", 3), 3, "separator");
    must(gzwrite(gz, qual, qual_len), qual_len, "quality");
    must(gzwrite(gz, "\n", 1), 1, "newline");
}

void FastqWriter::write_record_plain(std::FILE* fp, const std::string& sample_id, long long n,
                                     char suffix, const char* seq, int seq_len,
                                     const char* qual, int qual_len) {
    char header[256];
    int header_len = 0;
    if (suffix) {
        header_len = std::snprintf(header, sizeof(header), "@%s_%lld/%c\n",
                                   sample_id.c_str(), n, suffix);
    } else {
        header_len = std::snprintf(header, sizeof(header), "@%s_%lld\n",
                                   sample_id.c_str(), n);
    }
    if (header_len < 0 || header_len >= static_cast<int>(sizeof(header))) {
        throw std::runtime_error("FASTQ header overflow for sample " + sample_id);
    }

    auto must = [&](size_t written, size_t expected, const char* label) {
        if (written != expected) {
            throw std::runtime_error(std::string("fwrite failed for ") + label +
                                     " on sample " + sample_id + ": " + std::strerror(errno));
        }
    };

    must(std::fwrite(header, 1, header_len, fp), static_cast<size_t>(header_len), "header");
    must(std::fwrite(seq, 1, seq_len, fp), static_cast<size_t>(seq_len), "sequence");
    must(std::fwrite("\n+\n", 1, 3, fp), 3, "separator");
    must(std::fwrite(qual, 1, qual_len, fp), static_cast<size_t>(qual_len), "quality");
    must(std::fwrite("\n", 1, 1, fp), 1, "newline");
}

void FastqWriter::append(const std::string& sample_id, int lane,
                         const char* r1_seq, const char* r1_qual, int r1_len,
                         const char* r2_seq, const char* r2_qual, int r2_len) {
    const bool paired = (r2_len > 0 && r2_seq != nullptr);
    LaneStream& s = get_or_open(sample_id, lane, paired);

    const long long n = s.read_count;
    const char r1_suffix = paired ? '1' : '\0';

    if (gzip_) {
        write_record_gz(s.r1_gz, sample_id, n, r1_suffix, r1_seq, r1_len, r1_qual, r1_len);
        if (paired) {
            write_record_gz(s.r2_gz, sample_id, n, '2', r2_seq, r2_len, r2_qual, r2_len);
        }
    } else {
        write_record_plain(s.r1_plain, sample_id, n, r1_suffix, r1_seq, r1_len, r1_qual, r1_len);
        if (paired) {
            write_record_plain(s.r2_plain, sample_id, n, '2', r2_seq, r2_len, r2_qual, r2_len);
        }
    }

    s.read_count++;
}

void FastqWriter::close() {
    if (closed_) {
        return;
    }
    closed_ = true;

    std::string error_msg;
    for (auto& [_, s] : streams_) {
        if (s.r1_gz) {
            int rc = gzclose(s.r1_gz);
            if (rc != Z_OK && error_msg.empty()) {
                error_msg = "gzclose failed for an R1 stream";
            }
            s.r1_gz = nullptr;
        }
        if (s.r2_gz) {
            int rc = gzclose(s.r2_gz);
            if (rc != Z_OK && error_msg.empty()) {
                error_msg = "gzclose failed for an R2 stream";
            }
            s.r2_gz = nullptr;
        }
        if (s.r1_plain) {
            if (std::fflush(s.r1_plain) != 0 && error_msg.empty()) {
                error_msg = std::string("fflush failed: ") + std::strerror(errno);
            }
            std::fclose(s.r1_plain);
            s.r1_plain = nullptr;
        }
        if (s.r2_plain) {
            if (std::fflush(s.r2_plain) != 0 && error_msg.empty()) {
                error_msg = std::string("fflush failed: ") + std::strerror(errno);
            }
            std::fclose(s.r2_plain);
            s.r2_plain = nullptr;
        }
    }
    streams_.clear();

    if (!error_msg.empty()) {
        throw std::runtime_error(error_msg);
    }
}
