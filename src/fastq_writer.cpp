#include "fastq_writer.h"
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include "read_types.h"

// Dummy implementation: writes each barcode's reads to a FASTQ file in the output folder
void write_fastq(const std::string& output_folder, const std::unordered_map<std::string, std::vector<Read>>& demuxed_data) {
    // TODO: Replace with actual FASTQ writing logic
    for (const auto& [barcode, reads] : demuxed_data) {
        std::string filename = output_folder + "/" + barcode + ".fastq";
        std::ofstream out(filename);
        if (!out) throw std::runtime_error("Failed to open output file: " + filename);
        for (size_t i = 0; i < reads.size(); ++i) {
            out << "@read" << i << "\n";
            out << reads[i].sequence << "\n";
            out << "+\n";
            out << reads[i].quality << "\n";
        }
    }
}
