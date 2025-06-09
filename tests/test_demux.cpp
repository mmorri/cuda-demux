#include <cassert>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include "../include/demux.h"
#include "../include/read_types.h"

void test_demux_basic() {
    std::cout << "[TEST] test_demux_basic... ";
    // Construct a simple set of reads and barcodes
    std::vector<Read> reads = {
        {"ACGT", "IIII"},
        {"TGCA", "IIII"},
        {"ACGT", "IIII"}
    };
    // Simulate a sample sheet with barcodes
    std::string samplesheet = "test_samplesheet.csv";
    // For test: create a fake load_barcodes for this context
    // (In real code, mock load_barcodes or use dependency injection)
    // Here, just check that demux runs without error (integration test)
    try {
        auto result = demux(reads, samplesheet);
        assert(result.size() >= 0); // Should not throw
        std::cout << "PASSED\n";
    } catch (const std::exception& ex) {
        std::cerr << "FAILED: " << ex.what() << std::endl;
        assert(false);
    }
}

int main() {
    test_demux_basic();
    // Add more tests here
    std::cout << "[ALL TESTS PASSED]\n";
    return 0;
}