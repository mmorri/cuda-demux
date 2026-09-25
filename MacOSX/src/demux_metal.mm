#include "../include/metal_compute.h"
#include "../../include/common.h"
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <sstream>

using namespace metal_compute;

class DemuxMetal {
private:
    std::unique_ptr<MetalCompute> compute;
    std::shared_ptr<MetalBuffer> sample_barcodes_buffer;
    size_t num_samples;
    size_t barcode_length;
    size_t max_mismatches;

public:
    DemuxMetal(size_t max_mm = 1) : max_mismatches(max_mm) {
        compute = std::make_unique<MetalCompute>();

        // Load shader source
        std::string shader_path = "shaders/demux_kernels.metal";
        std::ifstream shader_file(shader_path);
        if (shader_file) {
            std::string shader_source((std::istreambuf_iterator<char>(shader_file)),
                                     std::istreambuf_iterator<char>());
            compute->loadShaderSource(shader_source);
        }

        // Create pipelines
        compute->createComputePipeline("match_barcodes");
        compute->createComputePipeline("trim_adapters");
    }

    void loadSampleBarcodes(const std::vector<std::string>& barcodes) {
        num_samples = barcodes.size();
        if (num_samples == 0) return;

        barcode_length = barcodes[0].length();

        // Flatten barcodes into contiguous array
        std::vector<char> flat_barcodes(num_samples * barcode_length);
        for (size_t i = 0; i < num_samples; i++) {
            memcpy(&flat_barcodes[i * barcode_length], barcodes[i].data(), barcode_length);
        }

        sample_barcodes_buffer = compute->createBuffer(flat_barcodes.size(), flat_barcodes.data());
    }

    std::vector<int> matchBarcodes(const std::vector<std::string>& read_barcodes) {
        size_t num_reads = read_barcodes.size();
        std::vector<int> matches(num_reads, -1);

        if (num_reads == 0 || num_samples == 0) return matches;

        // Flatten read barcodes
        std::vector<char> flat_reads(num_reads * barcode_length);
        for (size_t i = 0; i < num_reads; i++) {
            memcpy(&flat_reads[i * barcode_length], read_barcodes[i].data(), barcode_length);
        }

        auto reads_buffer = compute->createBuffer(flat_reads.size(), flat_reads.data());
        auto matches_buffer = compute->createBuffer(num_reads * sizeof(int));
        auto distances_buffer = compute->createBuffer(num_reads * sizeof(int));

        // Initialize distances to max
        std::vector<int> init_distances(num_reads, max_mismatches + 1);
        memcpy(distances_buffer->contents(), init_distances.data(), init_distances.size() * sizeof(int));

        std::vector<uint32_t> constants = {
            static_cast<uint32_t>(num_reads),
            static_cast<uint32_t>(num_samples),
            static_cast<uint32_t>(barcode_length),
            static_cast<uint32_t>(max_mismatches)
        };

        // Process in batches of samples
        size_t samples_per_thread = 8;
        size_t num_batches = (num_samples + samples_per_thread - 1) / samples_per_thread;

        MTLSize gridSize = MTLSizeMake(num_reads, num_batches, 1);
        MTLSize threadGroupSize = MTLSizeMake(32, 8, 1);

        compute->dispatch("match_barcodes",
                         {reads_buffer, sample_barcodes_buffer, matches_buffer, distances_buffer},
                         constants,
                         gridSize,
                         threadGroupSize);

        memcpy(matches.data(), matches_buffer->contents(), num_reads * sizeof(int));
        return matches;
    }

    std::vector<uint32_t> findAdapters(const std::vector<std::string>& sequences,
                                       const std::string& adapter,
                                       size_t min_overlap = 10) {
        size_t num_reads = sequences.size();
        if (num_reads == 0) return {};

        size_t read_length = sequences[0].length();
        std::vector<uint32_t> trim_positions(num_reads, read_length);

        // Flatten sequences
        std::vector<char> flat_seqs(num_reads * read_length);
        for (size_t i = 0; i < num_reads; i++) {
            memcpy(&flat_seqs[i * read_length], sequences[i].data(), read_length);
        }

        auto seqs_buffer = compute->createBuffer(flat_seqs.size(), flat_seqs.data());
        auto adapter_buffer = compute->createBuffer(adapter.size(), adapter.data());
        auto trim_buffer = compute->createBuffer(num_reads * sizeof(uint32_t));

        std::vector<uint32_t> constants = {
            static_cast<uint32_t>(num_reads),
            static_cast<uint32_t>(read_length),
            static_cast<uint32_t>(adapter.length()),
            static_cast<uint32_t>(min_overlap)
        };

        MTLSize gridSize = MTLSizeMake(num_reads, 1, 1);
        MTLSize threadGroupSize = MTLSizeMake(256, 1, 1);

        compute->dispatch("trim_adapters",
                         {seqs_buffer, adapter_buffer, trim_buffer},
                         constants,
                         gridSize,
                         threadGroupSize);

        memcpy(trim_positions.data(), trim_buffer->contents(), num_reads * sizeof(uint32_t));
        return trim_positions;
    }
};

// Main demux class that coordinates everything
class MetalDemultiplexer {
private:
    std::unique_ptr<DemuxMetal> demux;
    std::unordered_map<int, std::string> index_to_sample;
    std::unordered_map<std::string, std::ofstream> output_files;

public:
    MetalDemultiplexer(const std::string& samplesheet_path, size_t max_mismatches = 1) {
        demux = std::make_unique<DemuxMetal>(max_mismatches);
        loadSampleSheet(samplesheet_path);
    }

    void loadSampleSheet(const std::string& path) {
        std::ifstream file(path);
        std::string line;
        std::vector<std::string> barcodes;
        bool in_data_section = false;

        while (std::getline(file, line)) {
            if (line.find("[Data]") != std::string::npos) {
                in_data_section = true;
                std::getline(file, line); // Skip header
                continue;
            }

            if (in_data_section && !line.empty()) {
                std::istringstream iss(line);
                std::string sample_id, sample_name, index1, index2;

                // Parse CSV line
                std::getline(iss, sample_id, ',');
                std::getline(iss, sample_name, ',');
                // Skip some fields...
                for (int i = 0; i < 5; i++) {
                    std::string dummy;
                    std::getline(iss, dummy, ',');
                }
                std::getline(iss, index1, ',');
                std::getline(iss, index2, ',');

                std::string barcode = index1;
                if (!index2.empty() && index2 != "NA") {
                    barcode += index2;
                }

                int idx = barcodes.size();
                barcodes.push_back(barcode);
                index_to_sample[idx] = sample_name;
            }
        }

        demux->loadSampleBarcodes(barcodes);
    }

    void processReads(const std::vector<Read>& reads, const std::string& output_dir) {
        // Extract barcodes from reads
        std::vector<std::string> read_barcodes;
        for (const auto& read : reads) {
            read_barcodes.push_back(read.barcode_seq);
        }

        // Match barcodes using Metal
        auto matches = demux->matchBarcodes(read_barcodes);

        // Write to appropriate files
        for (size_t i = 0; i < reads.size(); i++) {
            int sample_idx = matches[i];
            std::string filename;

            if (sample_idx >= 0 && index_to_sample.find(sample_idx) != index_to_sample.end()) {
                filename = output_dir + "/" + index_to_sample[sample_idx] + "_R1.fastq";
            } else {
                filename = output_dir + "/Undetermined_R1.fastq";
            }

            if (output_files.find(filename) == output_files.end()) {
                output_files[filename].open(filename);
            }

            // Write FASTQ entry
            output_files[filename] << "@" << reads[i].header << "\n"
                                  << reads[i].sequence << "\n"
                                  << "+\n"
                                  << reads[i].quality << "\n";
        }
    }

    ~MetalDemultiplexer() {
        for (auto& [name, file] : output_files) {
            file.close();
        }
    }
};