#include "../include/metal_compute.h"
#include "../../include/bcl_parser.h"
#include <iostream>
#include <chrono>
#include <cstring>

using namespace metal_compute;

class BCLParserMetal {
private:
    std::unique_ptr<MetalCompute> compute;
    size_t max_batch_size;
    size_t current_batch_size;

public:
    BCLParserMetal() {
        compute = std::make_unique<MetalCompute>();

        // Load shader source (in production, would load from compiled .metallib)
        std::string shader_path = "shaders/demux_kernels.metal";
        std::ifstream shader_file(shader_path);
        if (shader_file) {
            std::string shader_source((std::istreambuf_iterator<char>(shader_file)),
                                     std::istreambuf_iterator<char>());
            compute->loadShaderSource(shader_source);
        }

        // Create pipelines for each kernel
        compute->createComputePipeline("convert_bcl_to_bases");
        compute->createComputePipeline("filter_by_quality");
        compute->createComputePipeline("reverse_complement");

        // Calculate max batch size based on available memory
        size_t available_memory = compute->getDeviceMemory();
        max_batch_size = available_memory / (1024 * 1024 * 10); // Conservative estimate
    }

    void processBCLBatch(const std::vector<uint8_t>& bcl_data,
                        std::vector<char>& output_bases,
                        std::vector<uint8_t>& output_quality,
                        size_t num_clusters,
                        size_t num_cycles) {

        // Prepare buffers
        auto bcl_buffer = compute->createBuffer(bcl_data.size(), bcl_data.data());
        auto bases_buffer = compute->createBuffer(output_bases.size());
        auto quality_buffer = compute->createBuffer(output_quality.size());

        // Setup constants
        std::vector<uint32_t> constants = {
            static_cast<uint32_t>(num_cycles),
            static_cast<uint32_t>(num_clusters)
        };

        // Calculate thread groups
        MTLSize gridSize = MTLSizeMake(num_clusters, num_cycles, 1);
        MTLSize threadGroupSize = MTLSizeMake(32, 32, 1); // Typical for Metal

        // Dispatch kernel
        compute->dispatch("convert_bcl_to_bases",
                         {bcl_buffer, bases_buffer, quality_buffer},
                         constants,
                         gridSize,
                         threadGroupSize);

        // Copy results back
        memcpy(output_bases.data(), bases_buffer->contents(), output_bases.size());
        memcpy(output_quality.data(), quality_buffer->contents(), output_quality.size());
    }

    void filterByQuality(const std::vector<uint8_t>& quality_scores,
                        std::vector<bool>& pass_filter,
                        size_t num_reads,
                        size_t read_length,
                        uint32_t min_quality,
                        float min_quality_fraction) {

        auto quality_buffer = compute->createBuffer(quality_scores.size(), quality_scores.data());
        auto filter_buffer = compute->createBuffer(pass_filter.size() * sizeof(bool));

        std::vector<uint32_t> constants = {
            static_cast<uint32_t>(num_reads),
            static_cast<uint32_t>(read_length),
            min_quality
        };

        // Note: float constant needs special handling
        float fraction = min_quality_fraction;
        constants.push_back(*reinterpret_cast<uint32_t*>(&fraction));

        MTLSize gridSize = MTLSizeMake(num_reads, 1, 1);
        MTLSize threadGroupSize = MTLSizeMake(256, 1, 1);

        compute->dispatch("filter_by_quality",
                         {quality_buffer, filter_buffer},
                         constants,
                         gridSize,
                         threadGroupSize);

        memcpy(pass_filter.data(), filter_buffer->contents(), pass_filter.size() * sizeof(bool));
    }
};

// C++ wrapper functions to match CUDA interface
extern "C" {

void* bcl_parser_metal_init() {
    return new BCLParserMetal();
}

void bcl_parser_metal_destroy(void* parser) {
    delete static_cast<BCLParserMetal*>(parser);
}

void bcl_parser_metal_process(void* parser,
                              const uint8_t* bcl_data,
                              char* output_bases,
                              uint8_t* output_quality,
                              size_t num_clusters,
                              size_t num_cycles) {
    auto p = static_cast<BCLParserMetal*>(parser);

    std::vector<uint8_t> bcl_vec(bcl_data, bcl_data + num_clusters * num_cycles);
    std::vector<char> bases_vec(num_clusters * num_cycles);
    std::vector<uint8_t> quality_vec(num_clusters * num_cycles);

    p->processBCLBatch(bcl_vec, bases_vec, quality_vec, num_clusters, num_cycles);

    memcpy(output_bases, bases_vec.data(), bases_vec.size());
    memcpy(output_quality, quality_vec.data(), quality_vec.size());
}

} // extern "C"