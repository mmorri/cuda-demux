#ifndef DEMUX_H
#define DEMUX_H

#include <string>
#include <vector>
#include "common.h"
#include "fastq_writer.h"

void demux_and_write(const std::vector<LaneBclData>& lanes,
                     const std::string& samplesheet,
                     const std::string& run_folder,
                     FastqWriter& writer);

std::vector<SampleInfo> load_sample_info(const std::string& samplesheet);
bool validate_sample_barcodes(const std::vector<SampleInfo>& samples);

#endif
