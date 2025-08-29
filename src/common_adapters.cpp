#include "common.h"

// Common Illumina adapter sequences that should be filtered out
// rather than treated as valid barcodes
const std::vector<std::string> COMMON_ADAPTERS = {
    "CTGTCTCTTATACACATCT", // Illumina Universal Adapter
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", // TruSeq Adapter, Read 1
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", // TruSeq Adapter, Read 2
    "AATGATACGGCGACCACCGAGATCTACAC",     // Illumina P5
    "CAAGCAGAAGACGGCATACGAGAT"           // Illumina P7
};
