#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <string>
#include <vector>

struct Read {
    std::string sequence;
    std::string quality;
    std::string index1;
    std::string index2;
    std::string read2_sequence;
    std::string read2_quality;
    int lane = 1;
};

struct SampleInfo {
    std::string sample_id;
    std::string index1;
    std::string index2;
    bool reverse_complement_i2 = false;
    int lane = 0;

    std::string getCombinedBarcode() const {
        if (reverse_complement_i2 && !index2.empty()) {
            return index1 + reverseComplement(index2);
        }
        return index1 + index2;
    }

    static std::string reverseComplement(const std::string& seq) {
        std::string rc;
        rc.reserve(seq.length());
        for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
            switch (*it) {
                case 'A': case 'a': rc += 'T'; break;
                case 'T': case 't': rc += 'A'; break;
                case 'G': case 'g': rc += 'C'; break;
                case 'C': case 'c': rc += 'G'; break;
                default: rc += 'N'; break;
            }
        }
        return rc;
    }
};

struct LaneBclData {
    int lane = 1;
    size_t num_clusters = 0;
    int total_cycles = 0;
    int r1_len = 0;
    int i1_len = 0;
    int i2_len = 0;
    int r2_len = 0;
    std::vector<int> read_segments;
    std::vector<std::vector<uint8_t>> bcl;
};

extern const std::vector<std::string> COMMON_ADAPTERS;

#endif
