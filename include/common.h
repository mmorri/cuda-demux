#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>

struct Read {
    std::string sequence;
    std::string quality;
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
