#include <iostream>
#include <string>
#include "demux.h"
#include "bcl_parser.h"
#include "fastq_writer.h"

int main(int argc, char* argv[]) {
    if (argc != 7) {
        std::cerr << "Usage: ./cuda-demux --input <BCL_FOLDER> --samplesheet <CSV> --output <OUTPUT_FOLDER>\n";
        return 1;
    }

    std::string input_folder = argv[2];
    std::string samplesheet = argv[4];
    std::string output_folder = argv[6];

    std::cout << "Parsing BCL files...\n";
    auto reads = parse_bcl(input_folder);

    std::cout << "Demultiplexing reads...\n";
    auto demuxed_data = demux(reads, samplesheet);

    std::cout << "Writing FASTQ files...\n";
    write_fastq(output_folder, demuxed_data);

    std::cout << "Demultiplexing completed successfully.\n";
    return 0;
}