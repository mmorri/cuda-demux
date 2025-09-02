#include <iostream>
#include <string>
#include <cstdlib>
#include "demux.h"
#include "bcl_parser.h"
#include "fastq_writer.h"

int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: ./cuda-demux --input <RUN_FOLDER> --samplesheet <CSV> --output <OUTPUT_FOLDER> [--gzip]\n";
        return 1;
    }

    std::string input_folder;
    std::string samplesheet;
    std::string output_folder;
    bool gzip_output = false;
    // Optional tuning parameters via CLI flags (mapped to env vars for simplicity)
    std::string opt_batch_size;
    std::string opt_mem_fraction;
    std::string opt_device;
    bool opt_no_adaptive = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--input" && i + 1 < argc) { input_folder = argv[++i]; }
        else if (arg == "--samplesheet" && i + 1 < argc) { samplesheet = argv[++i]; }
        else if (arg == "--output" && i + 1 < argc) { output_folder = argv[++i]; }
        else if (arg == "--gzip") { gzip_output = true; }
        else if (arg == "--batch-size" && i + 1 < argc) { opt_batch_size = argv[++i]; }
        else if (arg == "--gpu-mem-fraction" && i + 1 < argc) { opt_mem_fraction = argv[++i]; }
        else if (arg == "--device" && i + 1 < argc) { opt_device = argv[++i]; }
        else if (arg == "--no-adaptive-probe") { opt_no_adaptive = true; }
    }
    if (input_folder.empty() || samplesheet.empty() || output_folder.empty()) {
        std::cerr << "Error: Missing required arguments.\n";
        std::cerr << "Usage: ./cuda-demux --input <RUN_FOLDER> --samplesheet <CSV> --output <OUTPUT_FOLDER> [--gzip]" << std::endl;
        return 1;
    }

    // Map options to environment for downstream modules
    if (!opt_batch_size.empty()) setenv("CUDA_DEMUX_BATCH_SIZE", opt_batch_size.c_str(), 1);
    if (!opt_mem_fraction.empty()) setenv("CUDA_DEMUX_MEM_FRACTION", opt_mem_fraction.c_str(), 1);
    if (!opt_device.empty()) setenv("CUDA_DEMUX_DEVICE", opt_device.c_str(), 1);
    if (opt_no_adaptive) setenv("CUDA_DEMUX_NO_ADAPTIVE", "1", 1);

    std::cout << "Parsing BCL files...\n";
    auto reads = parse_bcl(input_folder);

    std::cout << "Demultiplexing reads...\n";
    auto demuxed_data = demux(reads, samplesheet, input_folder);

    std::cout << "Writing FASTQ files...\n";
    write_fastq(output_folder, demuxed_data, gzip_output);

    std::cout << "Demultiplexing completed successfully.\n";
    return 0;
}
