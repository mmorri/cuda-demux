#ifndef FASTQ_WRITER_H
#define FASTQ_WRITER_H

#include <cstdio>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <zlib.h>

class FastqWriter {
public:
    FastqWriter(std::string output_folder, bool gzip_output);
    ~FastqWriter();

    FastqWriter(const FastqWriter&) = delete;
    FastqWriter& operator=(const FastqWriter&) = delete;

    void append(const std::string& sample_id, int lane,
                const char* r1_seq, const char* r1_qual, int r1_len,
                const char* r2_seq, const char* r2_qual, int r2_len);

    void close();

    size_t samples_written() const { return distinct_samples_.size(); }

private:
    struct LaneStream {
        gzFile r1_gz = nullptr;
        gzFile r2_gz = nullptr;
        std::FILE* r1_plain = nullptr;
        std::FILE* r2_plain = nullptr;
        long long read_count = 0;
        bool paired = false;
    };

    LaneStream& get_or_open(const std::string& sample_id, int lane, bool paired);
    void write_record_gz(gzFile gz, const std::string& sample_id, long long n,
                         char suffix, const char* seq, int seq_len,
                         const char* qual, int qual_len);
    void write_record_plain(std::FILE* fp, const std::string& sample_id, long long n,
                            char suffix, const char* seq, int seq_len,
                            const char* qual, int qual_len);

    std::string folder_;
    bool gzip_;
    bool closed_ = false;
    std::unordered_map<std::string, LaneStream> streams_;
    std::unordered_set<std::string> distinct_samples_;
};

#endif
