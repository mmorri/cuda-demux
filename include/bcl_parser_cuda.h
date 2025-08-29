#ifndef BCL_PARSER_CUDA_H
#define BCL_PARSER_CUDA_H

#include "common.h"
#include <vector>

// This function is the public interface for the CUDA BCL decoding.
// It takes raw BCL data for a full tile across all cycles and decodes it on the GPU.
void decode_bcl_data_cuda(
    const std::vector<char*>& h_bcl_data, // Host pointers to raw BCL data for each cycle
    const std::vector<size_t>& h_bcl_sizes, // Sizes of each BCL buffer
    const std::vector<int>& h_read_structure, // Defines which cycle belongs to which read segment
    std::vector<Read>& reads, // Output vector of assembled reads
    size_t num_clusters
);

#endif // BCL_PARSER_CUDA_H
