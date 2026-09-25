#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>

#include "common.h"
#include "fastq_writer.h"

namespace fs = std::filesystem;

// CHECK() is compiled out in Release builds; these checks must always run.
#define CHECK(cond)                                                                   \
    do {                                                                              \
        if (!(cond)) {                                                                \
            std::cerr << __FILE__ << ":" << __LINE__ << ": check failed: " #cond "\n"; \
            std::exit(1);                                                             \
        }                                                                             \
    } while (0)

static fs::path make_temp_dir() {
    auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    fs::path dir = fs::temp_directory_path() / ("cuda-demux-test-" + std::to_string(stamp));
    fs::create_directories(dir);
    return dir;
}

static std::string read_text_file(const fs::path& path) {
    std::ifstream input(path);
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

static std::string read_gzip_file(const fs::path& path) {
    gzFile file = gzopen(path.string().c_str(), "rb");
    CHECK(file != nullptr);

    std::string out;
    char buffer[256];
    int bytes = 0;
    while ((bytes = gzread(file, buffer, sizeof(buffer))) > 0) {
        out.append(buffer, bytes);
    }
    gzclose(file);
    return out;
}

static void test_reverse_complement() {
    CHECK(SampleInfo::reverseComplement("ATGCN") == "NGCAT");
    CHECK(SampleInfo::reverseComplement("atgc") == "GCAT");
}

static void test_writer_single_end_plain() {
    fs::path out_dir = make_temp_dir();

    {
        FastqWriter w(out_dir.string(), false);
        const std::string seq1 = "ACGT";
        const std::string qual1 = "IIII";
        w.append("SampleA", 1, seq1.data(), qual1.data(), 4, nullptr, nullptr, 0);

        const std::string seq2 = "TGCA";
        const std::string qual2 = "HHHH";
        w.append("SampleA", 2, seq2.data(), qual2.data(), 4, nullptr, nullptr, 0);
        w.close();
    }

    CHECK(read_text_file(out_dir / "SampleA_L001_R1_001.fastq") ==
           "@SampleA_0\nACGT\n+\nIIII\n");
    CHECK(read_text_file(out_dir / "SampleA_L002_R1_001.fastq") ==
           "@SampleA_0\nTGCA\n+\nHHHH\n");

    fs::remove_all(out_dir);
}

static void test_writer_paired_end_gzip() {
    fs::path out_dir = make_temp_dir();

    {
        FastqWriter w(out_dir.string(), true);
        const std::string r1 = "AAAA";
        const std::string q1 = "IIII";
        const std::string r2 = "TTTT";
        const std::string q2 = "JJJJ";
        w.append("SampleB", 3, r1.data(), q1.data(), 4, r2.data(), q2.data(), 4);
        w.close();
    }

    CHECK(read_gzip_file(out_dir / "SampleB_L003_R1_001.fastq.gz") ==
           "@SampleB_0/1\nAAAA\n+\nIIII\n");
    CHECK(read_gzip_file(out_dir / "SampleB_L003_R2_001.fastq.gz") ==
           "@SampleB_0/2\nTTTT\n+\nJJJJ\n");

    fs::remove_all(out_dir);
}

static void test_writer_multiple_records() {
    fs::path out_dir = make_temp_dir();

    {
        FastqWriter w(out_dir.string(), false);
        const std::string seq1 = "ACGT";
        const std::string qual1 = "IIII";
        const std::string seq2 = "TTTT";
        const std::string qual2 = "JJJJ";
        w.append("SampleC", 1, seq1.data(), qual1.data(), 4, nullptr, nullptr, 0);
        w.append("SampleC", 1, seq2.data(), qual2.data(), 4, nullptr, nullptr, 0);
        w.close();
    }

    const std::string expected =
        "@SampleC_0\nACGT\n+\nIIII\n"
        "@SampleC_1\nTTTT\n+\nJJJJ\n";
    CHECK(read_text_file(out_dir / "SampleC_L001_R1_001.fastq") == expected);

    fs::remove_all(out_dir);
}

// Two samples whose (id, lane) pairs would collide under naive key concatenation.
static void test_writer_stream_key_collision() {
    fs::path out_dir = make_temp_dir();
    {
        FastqWriter w(out_dir.string(), false);
        w.append("Sample1", 1, "ACGT", "IIII", 4, nullptr, nullptr, 0);
        w.append("Sample", 11, "TTTT", "JJJJ", 4, nullptr, nullptr, 0);
        w.close();
    }
    CHECK(read_text_file(out_dir / "Sample1_L001_R1_001.fastq") == "@Sample1_0\nACGT\n+\nIIII\n");
    CHECK(read_text_file(out_dir / "Sample_L011_R1_001.fastq") == "@Sample_0\nTTTT\n+\nJJJJ\n");
    fs::remove_all(out_dir);
}

// A batch spanning several samples and many chunks: reads must land in the
// right file, in cluster order, with continuous numbering across batches, and
// the multi-member gzip output must decompress as one stream.
static void test_writer_batch() {
    fs::path out_dir = make_temp_dir();
    const int stride = 6;            // R1(2) + index(2) + R2(2)
    const size_t n = 70000;          // > 2 chunks of 32768
    std::vector<std::string> ids = {"A", "B", "undetermined"};
    std::string seq(n * stride, 'x'), qual(n * stride, 'q');
    std::vector<int> slot(n);
    for (size_t i = 0; i < n; ++i) {
        slot[i] = static_cast<int>(i % 3);
        seq[i * stride] = "ACGT"[i % 4];
        seq[i * stride + 4] = "TGCA"[i % 4];
        qual[i * stride] = static_cast<char>('!' + (i % 40));
    }
    FastqBatch b;
    b.lane = 2; b.count = n; b.seq = seq.data(); b.qual = qual.data(); b.stride = stride;
    b.r1_offset = 0; b.r1_len = 2; b.r2_offset = 4; b.r2_len = 2;
    b.sample_slot = slot.data(); b.sample_ids = &ids;
    {
        FastqWriter w(out_dir.string(), true);
        w.ensure_open("C", 2, true);     // never receives reads -> must be a valid empty gzip
        w.append_batch(b);
        w.append_batch(b);               // second batch continues the numbering
        w.close();
    }
    std::string a1 = read_gzip_file(out_dir / "A_L002_R1_001.fastq.gz");
    std::string a2 = read_gzip_file(out_dir / "A_L002_R2_001.fastq.gz");
    size_t records = 0;
    std::istringstream in(a1);
    std::string h, s, plus, q;
    size_t expect_i = 0;
    while (std::getline(in, h) && std::getline(in, s) && std::getline(in, plus) && std::getline(in, q)) {
        size_t i = (expect_i * 3) % n;         // cluster index of this record within the batch
        CHECK(h == "@A_" + std::to_string(records) + "/1");
        CHECK(s == std::string(1, "ACGT"[i % 4]) + "x");
        CHECK(q == std::string(1, static_cast<char>('!' + (i % 40))) + "q");
        ++records; expect_i = records % (n / 3 + (n % 3 > 0));
    }
    CHECK(records == 2 * ((n + 2) / 3));
    CHECK(a2.rfind("@A_0/2\nTx\n+\nqq\n", 0) == 0);
    CHECK(read_gzip_file(out_dir / "C_L002_R1_001.fastq.gz").empty());
    CHECK(fs::file_size(out_dir / "C_L002_R1_001.fastq.gz") > 0);
    fs::remove_all(out_dir);
}

int main() {
    test_reverse_complement();
    test_writer_single_end_plain();
    test_writer_paired_end_gzip();
    test_writer_multiple_records();
    test_writer_stream_key_collision();
    test_writer_batch();
    std::cout << "All unit tests passed.\n";
    return 0;
}
