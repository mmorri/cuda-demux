#ifndef DEMUX_H
#define DEMUX_H

#include <vector>
#include <unordered_map>
#include <string>

struct Read {
    std::string sequence;
    std::string quality;
};

std::unordered_map<std::string, std::vector<Read>> demux(const std::vector<Read>& reads, const std::string& samplesheet);

#endif // DEMUX_H