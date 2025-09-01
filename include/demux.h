#ifndef DEMUX_H
#define DEMUX_H

#include <vector>
#include <unordered_map>
#include "common.h"

// Demultiplex reads by sample, respecting SampleSheet lane restrictions and
// platform-specific i5 reverse-complement from RunParameters.xml.
// run_folder: path to Illumina run folder (to read RunParameters.xml)
std::unordered_map<std::string, std::vector<Read>> demux(
    const std::vector<Read>& reads,
    const std::string& samplesheet,
    const std::string& run_folder
);

#endif // DEMUX_H
