#include "bcl_parser.h"
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include "read_types.h"

// Dummy implementation: returns a vector of Read with dummy content for build success.
std::vector<Read> parse_bcl(const std::string& folder) {
    // TODO: Replace with actual BCL parsing logic
    std::vector<Read> reads;
    // Example dummy read
    reads.push_back(Read{"ACGTACGT", "IIIIIIII"});
    return reads;
}
