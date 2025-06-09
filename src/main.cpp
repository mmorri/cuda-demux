#include <iostream>
#include <string>
#include <exception>
#include "demux.h"
#include "bcl_parser.h"
#include "fastq_writer.h"
#include "read_types.h"

/**
 * @brief Entry point for CUDA-Demux.
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return Exit code.
 */
int main(int argc, char* argv[]) {
    try {
        if (argc != 7 || std::string(argv[1]) != "--input" || std::string(argv[3]) != "--samplesheet" || std::string(argv[5]) != "--output") {
            std::cerr << "Usage: ./cuda-demux --input <BCL_FOLDER> --samplesheet <CSV> --output <OUTPUT_FOLDER>\n";
            return 1;
        }

        const std::string input_folder = argv[2];
        const std::string samplesheet = argv[4];
        const std::string output_folder = argv[6];

        std::cout << "[INFO] Parsing BCL files...\n";
        auto reads = parse_bcl(input_folder);

        std::cout << "[INFO] Demultiplexing reads...\n";
        auto demuxed_data = demux(reads, samplesheet);

        std::cout << "[INFO] Writing FASTQ files...\n";
        write_fastq(output_folder, demuxed_data);

        std::cout << "[SUCCESS] Demultiplexing completed successfully.\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "[ERROR] " << ex.what() << std::endl;
        return 2;
    } catch (...) {
        std::cerr << "[ERROR] Unknown error occurred." << std::endl;
        return 3;
    }
}