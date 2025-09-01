#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>

struct Read {
    std::string sequence;     // Read 1 sequence
    std::string quality;      // Read 1 quality scores
    std::string index1;       // Index 1 sequence
    std::string index2;       // Index 2 sequence
    std::string read2_sequence;  // For paired-end reads
    std::string read2_quality;   // For paired-end reads
    int lane = 1;                // Lane number (1-8)
};

struct SampleInfo {
    std::string sample_id;
    std::string index1;
    std::string index2;
    bool reverse_complement_i2 = false;  // Platform-specific setting
    int lane = 0;                         // 0 means all lanes; 1-8 is specific lane
    
    // Helper to get the combined barcode
    std::string getCombinedBarcode() const {
        if (reverse_complement_i2 && !index2.empty()) {
            return index1 + reverseComplement(index2);
        }
        return index1 + index2;
    }
    
    // Helper to reverse complement a sequence
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

// Platform information from RunParameters.xml
struct PlatformInfo {
    std::string instrument_type;
    bool reverse_complement_i2 = false;
    int num_lanes = 1;
};

// List of common adapter sequences to filter out
extern const std::vector<std::string> COMMON_ADAPTERS;

#endif // COMMON_H
