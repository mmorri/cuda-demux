#ifndef GPU_DECOMPRESSOR_H
#define GPU_DECOMPRESSOR_H

#include <vector>
#include <cstddef>

// GPU-accelerated decompression functions
extern "C" {
    // Decompress a single block on GPU
    std::vector<char> decompress_block_gpu(
        const std::vector<char>& compressed_data,
        size_t expected_size
    );

    // Batch decompression for multiple blocks
    void decompress_batch_gpu(
        const std::vector<std::vector<char>>& compressed_blocks,
        std::vector<std::vector<char>>& decompressed_blocks,
        const std::vector<size_t>& expected_sizes
    );
}

#endif // GPU_DECOMPRESSOR_H