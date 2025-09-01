#include "fastq_writer.h"
#include "demux.h"
#include <fstream>
#include <iostream>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

namespace fs = std::filesystem;

void write_fastq(const std::string& output_folder, const std::unordered_map<std::string, std::vector<Read>>& demuxed_data, bool gzip_output) {
    try {
        fs::create_directories(output_folder);

        for (const auto& [sample_id, reads] : demuxed_data) {
            if (reads.empty()) continue;

            // Check if the data is paired-end by inspecting the first read
            bool is_paired_end = !reads.front().read2_sequence.empty();

            // Group reads by lane for lane-aware filenames
            std::unordered_map<int, std::vector<const Read*>> reads_by_lane;
            for (const auto& r : reads) {
                reads_by_lane[r.lane].push_back(&r);
            }

            for (const auto& [lane, lane_reads] : reads_by_lane) {
                char lane_buf[16];
                std::snprintf(lane_buf, sizeof(lane_buf), "L%03d", lane);
                std::string lane_tag(lane_buf);

                if (is_paired_end) {
                    std::string ext = gzip_output ? ".fastq.gz" : ".fastq";
                    std::string r1_filename = output_folder + "/" + sample_id + "_" + lane_tag + "_R1_001" + ext;
                    std::string r2_filename = output_folder + "/" + sample_id + "_" + lane_tag + "_R2_001" + ext;

                    if (gzip_output) {
                        gzFile r1_gz = gzopen(r1_filename.c_str(), "wb");
                        gzFile r2_gz = gzopen(r2_filename.c_str(), "wb");
                        if (!r1_gz || !r2_gz) {
                            std::cerr << "Error: Could not create gzip output files for sample " << sample_id << std::endl;
                            if (r1_gz) gzclose(r1_gz);
                            if (r2_gz) gzclose(r2_gz);
                            continue;
                        }
                        int read_count = 0;
                        for (const Read* rp : lane_reads) {
                            const Read& read = *rp;
                            std::string header1 = "@" + sample_id + "_" + std::to_string(read_count) + "/1\n";
                            std::string header2 = "@" + sample_id + "_" + std::to_string(read_count) + "/2\n";
                            gzwrite(r1_gz, header1.data(), header1.size());
                            gzwrite(r1_gz, read.sequence.data(), read.sequence.size()); gzwrite(r1_gz, "\n", 1);
                            gzwrite(r1_gz, "+\n", 2);
                            gzwrite(r1_gz, read.quality.data(), read.quality.size()); gzwrite(r1_gz, "\n", 1);

                            gzwrite(r2_gz, header2.data(), header2.size());
                            gzwrite(r2_gz, read.read2_sequence.data(), read.read2_sequence.size()); gzwrite(r2_gz, "\n", 1);
                            gzwrite(r2_gz, "+\n", 2);
                            gzwrite(r2_gz, read.read2_quality.data(), read.read2_quality.size()); gzwrite(r2_gz, "\n", 1);
                            read_count++;
                        }
                        gzclose(r1_gz);
                        gzclose(r2_gz);
                        std::cout << "Wrote " << lane_reads.size() << " paired-end reads to " << r1_filename << " and " << r2_filename << std::endl;
                    } else {
                        std::ofstream r1_outfile(r1_filename);
                        std::ofstream r2_outfile(r2_filename);
                        if (!r1_outfile.is_open() || !r2_outfile.is_open()) {
                            std::cerr << "Error: Could not create output files for sample " << sample_id << std::endl;
                            continue;
                        }
                        int read_count = 0;
                        for (const Read* rp : lane_reads) {
                            const Read& read = *rp;
                            r1_outfile << "@" << sample_id << "_" << read_count << "/1\n";
                            r1_outfile << read.sequence << "\n";
                            r1_outfile << "+\n";
                            r1_outfile << read.quality << "\n";
                            r2_outfile << "@" << sample_id << "_" << read_count << "/2\n";
                            r2_outfile << read.read2_sequence << "\n";
                            r2_outfile << "+\n";
                            r2_outfile << read.read2_quality << "\n";
                            read_count++;
                        }
                        r1_outfile.close();
                        r2_outfile.close();
                        std::cout << "Wrote " << lane_reads.size() << " paired-end reads to " << r1_filename << " and " << r2_filename << std::endl;
                    }
                } else {
                    std::string ext = gzip_output ? ".fastq.gz" : ".fastq";
                    std::string filename = output_folder + "/" + sample_id + "_" + lane_tag + "_R1_001" + ext;
                    if (gzip_output) {
                        gzFile gz = gzopen(filename.c_str(), "wb");
                        if (!gz) { std::cerr << "Error: Could not create gzip output file: " << filename << std::endl; continue; }
                        int read_count = 0;
                        for (const Read* rp : lane_reads) {
                            const Read& read = *rp;
                            std::string header = "@" + sample_id + "_" + std::to_string(read_count) + "\n";
                            gzwrite(gz, header.data(), header.size());
                            gzwrite(gz, read.sequence.data(), read.sequence.size()); gzwrite(gz, "\n", 1);
                            gzwrite(gz, "+\n", 2);
                            gzwrite(gz, read.quality.data(), read.quality.size()); gzwrite(gz, "\n", 1);
                            read_count++;
                        }
                        gzclose(gz);
                        std::cout << "Wrote " << lane_reads.size() << " single-end reads to " << filename << std::endl;
                    } else {
                        std::ofstream outfile(filename);
                        if (!outfile.is_open()) { std::cerr << "Error: Could not create output file: " << filename << std::endl; continue; }
                        int read_count = 0;
                        for (const Read* rp : lane_reads) {
                            const Read& read = *rp;
                            outfile << "@" << sample_id << "_" << read_count << "\n";
                            outfile << read.sequence << "\n";
                            outfile << "+\n";
                            outfile << read.quality << "\n";
                            read_count++;
                        }
                        outfile.close();
                        std::cout << "Wrote " << lane_reads.size() << " single-end reads to " << filename << std::endl;
                    }
                }
            }
        }
        std::cout << "Successfully wrote FASTQ files for " << demuxed_data.size() << " samples." << std::endl;
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error while writing FASTQ files: " << e.what() << std::endl;
    }
}
