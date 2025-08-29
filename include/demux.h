#ifndef DEMUX_H
#define DEMUX_H

#include <vector>
#include <unordered_map>
#include "common.h"

std::unordered_map<std::string, std::vector<Read>> demux(const std::vector<Read>& reads, const std::string& samplesheet);

#endif // DEMUX_H