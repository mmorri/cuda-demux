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
    // Create output directory if it doesn't exist
    try {
        if (!fs::exists(output_folder)) {
            fs::create_directories(output_folder);
        }
        
        // Write each sample's reads to a separate FASTQ file
        for (const auto& [barcode, reads] : demuxed_data) {
            // Create a filename based on the barcode
            std::string filename = output_folder + "/" + barcode + ".fastq";
            std::ofstream outfile(filename);
            
            if (!outfile.is_open()) {
                std::cerr << "Error: Could not create output file: " << filename << std::endl;
                continue;
            }
            
            // Write each read to the FASTQ file
            int read_count = 0;
            for (const auto& read : reads) {
                // FASTQ format:
                // Line 1: @read_id
                // Line 2: sequence
                // Line 3: +
                // Line 4: quality scores
                outfile << "@read_" << barcode << "_" << read_count << "\n";
                outfile << read.sequence << "\n";
                outfile << "+\n";
                outfile << read.quality << "\n";
                
                read_count++;
            }
            
            outfile.close();
            std::cout << "Wrote " << read_count << " reads to " << filename << std::endl;
        }
        
        std::cout << "Successfully wrote FASTQ files for " << demuxed_data.size() << " samples" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error writing FASTQ files: " << e.what() << std::endl;
    }
}
