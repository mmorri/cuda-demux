#include "fastq_writer.h"
#include "demux.h"
#include <fstream>
#include <iostream>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

void write_fastq(const std::string& output_folder, const std::unordered_map<std::string, std::vector<Read>>& demuxed_data) {
    try {
        fs::create_directories(output_folder);

        for (const auto& [sample_id, reads] : demuxed_data) {
            if (reads.empty()) continue;

            // Check if the data is paired-end by inspecting the first read
            bool is_paired_end = !reads.front().read2_sequence.empty();

            if (is_paired_end) {
                // Paired-end: write R1 and R2 files
                std::string r1_filename = output_folder + "/" + sample_id + "_R1.fastq";
                std::string r2_filename = output_folder + "/" + sample_id + "_R2.fastq";
                std::ofstream r1_outfile(r1_filename);
                std::ofstream r2_outfile(r2_filename);

                if (!r1_outfile.is_open() || !r2_outfile.is_open()) {
                    std::cerr << "Error: Could not create output files for sample " << sample_id << std::endl;
                    continue;
                }

                int read_count = 0;
                for (const auto& read : reads) {
                    // Write Read 1
                    r1_outfile << "@" << sample_id << "_" << read_count << "/1\n";
                    r1_outfile << read.sequence << "\n";
                    r1_outfile << "+\n";
                    r1_outfile << read.quality << "\n";

                    // Write Read 2
                    r2_outfile << "@" << sample_id << "_" << read_count << "/2\n";
                    r2_outfile << read.read2_sequence << "\n";
                    r2_outfile << "+\n";
                    r2_outfile << read.read2_quality << "\n";
                    
                    read_count++;
                }
                r1_outfile.close();
                r2_outfile.close();
                std::cout << "Wrote " << read_count << " paired-end reads to " << r1_filename << " and " << r2_filename << std::endl;

            } else {
                // Single-end: write one file
                std::string filename = output_folder + "/" + sample_id + ".fastq";
                std::ofstream outfile(filename);

                if (!outfile.is_open()) {
                    std::cerr << "Error: Could not create output file: " << filename << std::endl;
                    continue;
                }

                int read_count = 0;
                for (const auto& read : reads) {
                    outfile << "@" << sample_id << "_" << read_count << "\n";
                    outfile << read.sequence << "\n";
                    outfile << "+\n";
                    outfile << read.quality << "\n";
                    read_count++;
                }
                outfile.close();
                std::cout << "Wrote " << read_count << " single-end reads to " << filename << std::endl;
            }
        }
        std::cout << "Successfully wrote FASTQ files for " << demuxed_data.size() << " samples." << std::endl;
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error while writing FASTQ files: " << e.what() << std::endl;
    }
}
