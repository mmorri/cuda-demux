#ifndef DEMUX_H
#define DEMUX_H

#include <memory>
#include <string>
#include <vector>
#include "bcl_parser.h"
#include "common.h"
#include "fastq_writer.h"

// GPU demultiplexer for one run. Construct once (loads the sample sheet,
// selects the GPU, allocates the pipeline buffers), then feed lanes one at a
// time so only a single lane needs to be resident in host memory.
class Demuxer {
public:
    Demuxer(const RunLayout& run, const std::string& samplesheet, FastqWriter& writer);
    ~Demuxer();
    Demuxer(const Demuxer&) = delete;
    Demuxer& operator=(const Demuxer&) = delete;

    void process_lane(const LaneBclData& lane);
    void print_summary() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

std::vector<SampleInfo> load_sample_info(const std::string& samplesheet);
bool validate_sample_barcodes(const std::vector<SampleInfo>& samples);

#endif
