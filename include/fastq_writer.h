#ifndef FASTQ_WRITER_H
#define FASTQ_WRITER_H

#include <cstdio>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// A batch of decoded reads laid out as fixed-stride rows (one row per cluster)
// with R1 and R2 at fixed offsets inside each row.
struct FastqBatch {
    int lane = 1;
    size_t count = 0;
    const char* seq = nullptr;      // count * stride bytes
    const char* qual = nullptr;     // count * stride bytes
    size_t stride = 0;
    int r1_offset = 0;
    int r1_len = 0;
    int r2_offset = 0;
    int r2_len = 0;                 // 0 -> single-end
    const int* sample_slot = nullptr;                   // per read, index into *sample_ids
    const std::vector<std::string>* sample_ids = nullptr;
};

class FastqWriter {
public:
    FastqWriter(std::string output_folder, bool gzip_output, int gzip_level = 1);
    ~FastqWriter();

    FastqWriter(const FastqWriter&) = delete;
    FastqWriter& operator=(const FastqWriter&) = delete;

    // Creates the output file(s) for a sample/lane up front so every sample gets
    // a (possibly empty) FASTQ even when no read matches it.
    void ensure_open(const std::string& sample_id, int lane, bool paired);

    // Writes a whole batch. Reads are grouped by sample, formatted and (when
    // gzip is enabled) compressed in parallel, then appended in cluster order.
    // Each compressed chunk is an independent gzip member, so the output is a
    // standard multi-member gzip stream.
    void append_batch(const FastqBatch& batch);

    // Convenience single-read path (used by tests); equivalent to a 1-read batch.
    void append(const std::string& sample_id, int lane,
                const char* r1_seq, const char* r1_qual, int r1_len,
                const char* r2_seq, const char* r2_qual, int r2_len);

    void close();

    size_t samples_written() const { return distinct_samples_.size(); }

private:
    struct LaneStream {
        std::FILE* r1 = nullptr;
        std::FILE* r2 = nullptr;
        long long read_count = 0;
        bool paired = false;
    };

    LaneStream& get_or_open(const std::string& sample_id, int lane, bool paired);
    void write_all(std::FILE* fp, const std::string& data, const std::string& what);

    std::string folder_;
    bool gzip_;
    int gzip_level_;
    bool closed_ = false;
    std::unordered_map<std::string, LaneStream> streams_;
    std::unordered_set<std::string> distinct_samples_;
    // Reused across batches so large buffers are not re-faulted every batch.
    std::vector<std::string> chunk_pool_;     // one per task
    std::vector<std::string> format_pool_;    // one per thread
};

#endif
