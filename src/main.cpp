#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>

#include "bcl_parser.h"
#include "demux.h"
#include "fastq_writer.h"

int main(int argc, char* argv[]) {
    try {
        std::string input_folder;
        std::string samplesheet;
        std::string output_folder;
        bool gzip_output = false;
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
            std::cerr << "Usage: cuda-demux --input <RUN_FOLDER> --samplesheet <CSV> "
                         "--output <OUTPUT_FOLDER> [--gzip]\n"
                         "  [--batch-size N] [--gpu-mem-fraction F] [--device N] "
                         "[--no-adaptive-probe]" << std::endl;
            return 1;
        }
        if (!opt_batch_size.empty()) setenv("CUDA_DEMUX_BATCH_SIZE", opt_batch_size.c_str(), 1);
        if (!opt_mem_fraction.empty()) setenv("CUDA_DEMUX_MEM_FRACTION", opt_mem_fraction.c_str(), 1);
        if (!opt_device.empty()) setenv("CUDA_DEMUX_DEVICE", opt_device.c_str(), 1);
        if (opt_no_adaptive) setenv("CUDA_DEMUX_NO_ADAPTIVE", "1", 1);

        std::cout << "Parsing BCL files..." << std::endl;
        std::vector<LaneBclData> lanes = parse_bcl(input_folder);
        if (lanes.empty()) {
            std::cerr << "No lanes parsed; aborting." << std::endl;
            return 2;
        }

        FastqWriter writer(output_folder, gzip_output);
        std::cout << "Demultiplexing and writing FASTQ files..." << std::endl;
        demux_and_write(lanes, samplesheet, input_folder, writer);
        writer.close();

        std::cout << "Demultiplexing completed successfully." << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Fatal error: unknown exception" << std::endl;
        return 1;
    }
}
