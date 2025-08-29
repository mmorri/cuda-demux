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
};

struct SampleInfo {
    std::string sample_id;
    std::string index1;
    std::string index2;
    
    // Helper to get the combined barcode
    std::string getCombinedBarcode() const {
        return index1 + index2;
    }
};

// List of common adapter sequences to filter out
extern const std::vector<std::string> COMMON_ADAPTERS;

#endif // COMMON_H
