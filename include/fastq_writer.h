#ifndef FASTQ_WRITER_H
#define FASTQ_WRITER_H

#include <unordered_map>
#include <vector>
#include <string>
#include "read_types.h"

/**
 * @brief Write demultiplexed reads to FASTQ files.
 * @param output_folder Path to output directory.
 * @param demuxed_data Map from barcode/sample to reads.
 * @throws std::runtime_error on file errors.
 */
void write_fastq(const std::string& output_folder, const std::unordered_map<std::string, std::vector<Read>>& demuxed_data);

#endif // FASTQ_WRITER_H