#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "demux.h"

// Simple CSV barcode loader (assumes barcode is first column)
std::vector<std::string> load_barcodes(const std::string& samplesheet) {
    std::vector<std::string> barcodes;
    std::ifstream file(samplesheet);
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string barcode;
        if (std::getline(ss, barcode, ',')) {
            barcodes.push_back(barcode);
        }
    }
    return barcodes;
}
