#ifndef FASTQ_WRITER_H
#define FASTQ_WRITER_H

#include <unordered_map>
#include <vector>
#include "common.h"

void write_fastq(const std::string& output_folder, const std::unordered_map<std::string, std::vector<Read>>& demuxed_data);

#endif // FASTQ_WRITER_H