#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <zlib.h>

#include "common.h"
#include "fastq_writer.h"

namespace fs = std::filesystem;

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
    assert(file != nullptr);

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
    assert(SampleInfo::reverseComplement("ATGCN") == "NGCAT");
    assert(SampleInfo::reverseComplement("atgc") == "GCAT");
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

    assert(read_text_file(out_dir / "SampleA_L001_R1_001.fastq") ==
           "@SampleA_0\nACGT\n+\nIIII\n");
    assert(read_text_file(out_dir / "SampleA_L002_R1_001.fastq") ==
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

    assert(read_gzip_file(out_dir / "SampleB_L003_R1_001.fastq.gz") ==
           "@SampleB_0/1\nAAAA\n+\nIIII\n");
    assert(read_gzip_file(out_dir / "SampleB_L003_R2_001.fastq.gz") ==
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
    assert(read_text_file(out_dir / "SampleC_L001_R1_001.fastq") == expected);

    fs::remove_all(out_dir);
}

int main() {
    test_reverse_complement();
    test_writer_single_end_plain();
    test_writer_paired_end_gzip();
    test_writer_multiple_records();
    std::cout << "All unit tests passed.\n";
    return 0;
}
