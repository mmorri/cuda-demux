#ifndef DEMUX_H
#define DEMUX_H

#include <vector>
#include <unordered_map>
#include <string>
#include "read_types.h"

// Loads barcodes from a sample sheet CSV (assumes barcode is first column)
std::vector<std::string> load_barcodes(const std::string& samplesheet);

/**
 * @brief Demultiplex reads using a sample sheet.
 * @param reads Vector of sequencing reads.
 * @param samplesheet Path to the sample sheet CSV.
 * @return Map from barcode/sample to vector of reads.
 * @throws std::runtime_error on error.
 */
std::unordered_map<std::string, std::vector<Read>> demux(const std::vector<Read>& reads, const std::string& samplesheet);

#endif // DEMUX_H