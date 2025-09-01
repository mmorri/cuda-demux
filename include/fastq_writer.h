#ifndef FASTQ_WRITER_H
#define FASTQ_WRITER_H

#include <unordered_map>
#include <vector>
#include "common.h"

// gzip_output: if true, write .fastq.gz using zlib; else .fastq
void write_fastq(const std::string& output_folder,
                 const std::unordered_map<std::string, std::vector<Read>>& demuxed_data,
                 bool gzip_output);

#endif // FASTQ_WRITER_H
