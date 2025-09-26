#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <chrono>
#include "../../include/common.h"
#include "../../include/bcl_parser.h"
#include "../../include/fastq_writer.h"

namespace fs = std::filesystem;

// Forward declarations for Metal implementations
extern "C" {
    void* bcl_parser_metal_init();
    void bcl_parser_metal_destroy(void* parser);
    void bcl_parser_metal_process(void* parser,
                                  const uint8_t* bcl_data,
                                  char* output_bases,
                                  uint8_t* output_quality,
                                  size_t num_clusters,
                                  size_t num_cycles);
}

class MetalDemultiplexer;

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS]\n"
              << "Options:\n"
              << "  --input <path>         Input BCL/CBCL directory\n"
              << "  --output <path>        Output directory for FASTQ files\n"
              << "  --samplesheet <path>   Sample sheet CSV file\n"
              << "  --threads <n>          Number of CPU threads (default: auto)\n"
              << "  --batch-size <n>       GPU batch size (default: auto)\n"
              << "  --max-mismatches <n>   Maximum barcode mismatches (default: 1)\n"
              << "  --min-quality <n>      Minimum quality score (default: 20)\n"
              << "  --adapter <seq>        Adapter sequence to trim\n"
              << "  --help                 Show this help message\n";
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    std::string input_dir, output_dir, samplesheet;
    int num_threads = std::thread::hardware_concurrency();
    size_t batch_size = 0; // Auto
    int max_mismatches = 1;
    int min_quality = 20;
    std::string adapter_seq;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--input" && i + 1 < argc) {
            input_dir = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (arg == "--samplesheet" && i + 1 < argc) {
            samplesheet = argv[++i];
        } else if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (arg == "--batch-size" && i + 1 < argc) {
            batch_size = std::stoull(argv[++i]);
        } else if (arg == "--max-mismatches" && i + 1 < argc) {
            max_mismatches = std::stoi(argv[++i]);
        } else if (arg == "--min-quality" && i + 1 < argc) {
            min_quality = std::stoi(argv[++i]);
        } else if (arg == "--adapter" && i + 1 < argc) {
            adapter_seq = argv[++i];
        } else if (arg == "--help") {
            printUsage(argv[0]);
            return 0;
        }
    }

    // Validate required arguments
    if (input_dir.empty() || output_dir.empty() || samplesheet.empty()) {
        std::cerr << "Error: Missing required arguments\n\n";
        printUsage(argv[0]);
        return 1;
    }

    // Check Metal support
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) {
            std::cerr << "Error: Metal is not supported on this Mac\n";
            std::cerr << "This application requires a Mac with Metal support\n";
            return 1;
        }

        std::cout << "Using Metal device: " << [[device name] UTF8String] << "\n";
        std::cout << "Max threads per threadgroup: " << [device maxThreadsPerThreadgroup].width << "\n";
        std::cout << "Recommended working set: "
                  << ([device recommendedMaxWorkingSetSize] / (1024.0 * 1024.0 * 1024.0))
                  << " GB\n\n";
    }

    try {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Create output directory
        fs::create_directories(output_dir);

        // Initialize Metal BCL parser
        void* metal_parser = bcl_parser_metal_init();
        if (!metal_parser) {
            throw std::runtime_error("Failed to initialize Metal compute");
        }

        // Initialize CPU BCL parser for file I/O
        BCLParser bcl_parser(input_dir);

        // Parse RunInfo.xml
        RunInfo run_info = bcl_parser.parseRunInfo();
        std::cout << "Run Info:\n"
                  << "  Instrument: " << run_info.instrument << "\n"
                  << "  Run ID: " << run_info.run_id << "\n"
                  << "  Flow Cell: " << run_info.flowcell_id << "\n"
                  << "  Total Cycles: " << run_info.total_cycles << "\n"
                  << "  Lanes: " << run_info.lanes.size() << "\n\n";

        // Process tiles
        size_t total_reads = 0;
        size_t total_clusters = 0;

        for (const auto& lane : run_info.lanes) {
            std::cout << "Processing Lane " << lane << "...\n";

            for (const auto& surface : run_info.surfaces) {
                for (const auto& swath : run_info.swaths) {
                    for (const auto& tile : run_info.tiles) {
                        std::string tile_id = std::to_string(lane) + "_" +
                                            std::to_string(surface) +
                                            std::to_string(swath) +
                                            std::to_string(tile);

                        // Read BCL data for this tile
                        std::vector<std::vector<uint8_t>> cycles_data;
                        for (int cycle = 1; cycle <= run_info.total_cycles; cycle++) {
                            auto bcl_data = bcl_parser.readBCLFile(lane, tile_id, cycle);
                            cycles_data.push_back(bcl_data);
                        }

                        if (cycles_data.empty() || cycles_data[0].empty()) {
                            continue;
                        }

                        size_t num_clusters = cycles_data[0].size();
                        total_clusters += num_clusters;

                        // Prepare flattened data for Metal processing
                        std::vector<uint8_t> flat_bcl(num_clusters * run_info.total_cycles);
                        for (size_t cycle = 0; cycle < run_info.total_cycles; cycle++) {
                            for (size_t cluster = 0; cluster < num_clusters; cluster++) {
                                flat_bcl[cycle * num_clusters + cluster] = cycles_data[cycle][cluster];
                            }
                        }

                        // Process with Metal
                        std::vector<char> bases(num_clusters * run_info.total_cycles);
                        std::vector<uint8_t> qualities(num_clusters * run_info.total_cycles);

                        bcl_parser_metal_process(metal_parser,
                                               flat_bcl.data(),
                                               bases.data(),
                                               qualities.data(),
                                               num_clusters,
                                               run_info.total_cycles);

                        // Convert to reads
                        std::vector<Read> reads;
                        for (size_t i = 0; i < num_clusters; i++) {
                            Read read;
                            read.header = tile_id + ":" + std::to_string(i);
                            read.sequence.assign(&bases[i * run_info.total_cycles], run_info.total_cycles);
                            read.quality.assign(&qualities[i * run_info.total_cycles], run_info.total_cycles);

                            // Extract barcode (simplified - should use read structure from RunInfo)
                            if (run_info.index1_cycles > 0) {
                                size_t barcode_start = run_info.read1_cycles;
                                read.barcode_seq = read.sequence.substr(barcode_start, run_info.index1_cycles);
                            }

                            reads.push_back(read);
                        }

                        total_reads += reads.size();

                        // TODO: Demultiplex and write to files
                        // This would use the MetalDemultiplexer class from demux_metal.mm
                    }
                }
            }
        }

        // Cleanup
        bcl_parser_metal_destroy(metal_parser);

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

        std::cout << "\nProcessing Complete!\n"
                  << "Total reads processed: " << total_reads << "\n"
                  << "Total clusters: " << total_clusters << "\n"
                  << "Time elapsed: " << duration.count() << " seconds\n"
                  << "Throughput: " << (total_reads / (duration.count() + 1)) << " reads/sec\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}