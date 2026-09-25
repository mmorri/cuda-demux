#include <chrono>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

#include "bcl_parser.h"
#include "demux.h"
#include "fastq_writer.h"

namespace {

void print_usage() {
    std::cerr << "Usage: cuda-demux --input <RUN_FOLDER> --samplesheet <CSV> --output <OUTPUT_FOLDER>\n"
                 "  [--gzip] [--gzip-level N]        gzip FASTQ output (level 1-9, default 1)\n"
                 "  [--barcode-mismatches N]         mismatches allowed per index (0-4, default 1)\n"
                 "  [--batch-size N] [--gpu-mem-fraction F] [--device N]\n"
                 "  [--verbose]\n";
}

double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        std::string input_folder;
        std::string samplesheet;
        std::string output_folder;
        bool gzip_output = false;
        int gzip_level = 1;
        std::string opt_batch_size;
        std::string opt_mem_fraction;
        std::string opt_device;
        std::string opt_mismatches;

        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            const bool has_value = (i + 1 < argc);
            if (arg == "--input" && has_value) { input_folder = argv[++i]; }
            else if (arg == "--samplesheet" && has_value) { samplesheet = argv[++i]; }
            else if (arg == "--output" && has_value) { output_folder = argv[++i]; }
            else if (arg == "--gzip") { gzip_output = true; }
            else if (arg == "--gzip-level" && has_value) { gzip_output = true; gzip_level = std::stoi(argv[++i]); }
            else if (arg == "--barcode-mismatches" && has_value) { opt_mismatches = argv[++i]; }
            else if (arg == "--batch-size" && has_value) { opt_batch_size = argv[++i]; }
            else if (arg == "--gpu-mem-fraction" && has_value) { opt_mem_fraction = argv[++i]; }
            else if (arg == "--device" && has_value) { opt_device = argv[++i]; }
            else if (arg == "--no-adaptive-probe") { /* accepted for compatibility; no effect */ }
            else if (arg == "--verbose") { setenv("CUDA_DEMUX_VERBOSE", "1", 1); }
            else if (arg == "--help" || arg == "-h") { print_usage(); return 0; }
            else {
                std::cerr << "Unknown or incomplete argument: " << arg << "\n";
                print_usage();
                return 1;
            }
        }
        if (input_folder.empty() || samplesheet.empty() || output_folder.empty()) {
            print_usage();
            return 1;
        }
        if (!opt_batch_size.empty()) setenv("CUDA_DEMUX_BATCH_SIZE", opt_batch_size.c_str(), 1);
        if (!opt_mem_fraction.empty()) setenv("CUDA_DEMUX_MEM_FRACTION", opt_mem_fraction.c_str(), 1);
        if (!opt_device.empty()) setenv("CUDA_DEMUX_DEVICE", opt_device.c_str(), 1);
        if (!opt_mismatches.empty()) setenv("CUDA_DEMUX_MISMATCHES", opt_mismatches.c_str(), 1);

        const auto t_start = std::chrono::steady_clock::now();
        const RunLayout run = parse_run_layout(input_folder);

        FastqWriter writer(output_folder, gzip_output, gzip_level);
        Demuxer demuxer(run, samplesheet, writer);

        // One lane at a time: parse, demux, release.
        double t_parse = 0.0, t_demux = 0.0;
        size_t lanes_done = 0;
        for (size_t i = 0; i < run.lane_dirs.size(); ++i) {
            const auto t0 = std::chrono::steady_clock::now();
            const LaneBclData lane = parse_lane(run, i);
            const auto t1 = std::chrono::steady_clock::now();
            if (lane.num_clusters == 0) continue;
            demuxer.process_lane(lane);
            const auto t2 = std::chrono::steady_clock::now();
            t_parse += std::chrono::duration<double>(t1 - t0).count();
            t_demux += std::chrono::duration<double>(t2 - t1).count();
            ++lanes_done;
        }
        if (lanes_done == 0) {
            std::cerr << "No lanes with passing-filter clusters; aborting." << std::endl;
            return 2;
        }
        writer.close();
        demuxer.print_summary();
        std::cout << std::fixed << std::setprecision(1)
                  << "BCL parsing took " << t_parse << " s, demultiplexing took " << t_demux
                  << " s (total " << seconds_since(t_start) << " s)" << std::endl;

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
