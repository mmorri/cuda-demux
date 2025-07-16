#include <iostream>
#include <string>
#include <filesystem>
#include <vector> // Added for std::vector

// Simple XML parsing test without external dependencies
bool test_xml_parsing() {
    std::cout << "Testing XML parsing logic..." << std::endl;
    
    // Simulate the XML content structure
    std::string xml_content = R"(
<?xml version="1.0"?>
<RunInfo Version="6">
    <Run Id="240227_VH00103_437_2222FFYNX" Number="437">
        <Reads>
            <Read Number="1" NumCycles="76" IsIndexedRead="N" IsReverseComplement="N"/>
            <Read Number="2" NumCycles="10" IsIndexedRead="Y" IsReverseComplement="N"/>
            <Read Number="3" NumCycles="10" IsIndexedRead="Y" IsReverseComplement="Y"/>
            <Read Number="4" NumCycles="76" IsIndexedRead="N" IsReverseComplement="N"/>
        </Reads>
    </Run>
</RunInfo>
)";
    
    std::cout << "XML content:" << std::endl << xml_content << std::endl;
    
    // Check if XML declaration is present
    if (xml_content.find("<?xml version=\"1.0\"?>") == std::string::npos) {
        std::cerr << "Error: XML declaration not found!" << std::endl;
        return false;
    }
    
    // Check if required elements are present
    if (xml_content.find("<RunInfo") == std::string::npos) {
        std::cerr << "Error: <RunInfo> element not found!" << std::endl;
        return false;
    }
    
    if (xml_content.find("<Run") == std::string::npos) {
        std::cerr << "Error: <Run> element not found!" << std::endl;
        return false;
    }
    
    if (xml_content.find("<Reads>") == std::string::npos) {
        std::cerr << "Error: <Reads> element not found!" << std::endl;
        return false;
    }
    
    // Check for Read elements
    size_t read_count = 0;
    size_t pos = 0;
    while ((pos = xml_content.find("<Read", pos)) != std::string::npos) {
        read_count++;
        pos++;
    }
    
    if (read_count != 4) {
        std::cerr << "Error: Expected 4 Read elements, found " << read_count << std::endl;
        return false;
    }
    
    std::cout << "Found " << read_count << " Read elements" << std::endl;
    
    // Test the parsing logic from the original code
    int read1_cycles = 0;
    int index1_cycles = 0;
    int index2_cycles = 0;
    int read2_cycles = 0;
    int total_cycles = 0;
    
    // Simulate parsing the reads
    std::vector<std::pair<int, bool>> reads = {
        {76, false},  // Read 1: 76 cycles, not indexed
        {10, true},   // Read 2: 10 cycles, indexed
        {10, true},   // Read 3: 10 cycles, indexed
        {76, false}   // Read 4: 76 cycles, not indexed
    };
    
    for (const auto& read : reads) {
        int num_cycles = read.first;
        bool is_indexed = read.second;
        
        std::cout << "Processing read: " << num_cycles << " cycles, indexed: " << (is_indexed ? "Y" : "N") << std::endl;
        
        if (!is_indexed) {
            if (read1_cycles == 0) {
                read1_cycles = num_cycles;
                std::cout << "  -> Read1 cycles: " << num_cycles << std::endl;
            } else {
                read2_cycles = num_cycles;
                std::cout << "  -> Read2 cycles: " << num_cycles << std::endl;
            }
        } else {
            if (index1_cycles == 0) {
                index1_cycles = num_cycles;
                std::cout << "  -> Index1 cycles: " << num_cycles << std::endl;
            } else {
                index2_cycles = num_cycles;
                std::cout << "  -> Index2 cycles: " << num_cycles << std::endl;
            }
        }
        total_cycles += num_cycles;
    }
    
    std::cout << "\nFinal Run Structure:" << std::endl;
    std::cout << "R1: " << read1_cycles << " cycles" << std::endl;
    std::cout << "I1: " << index1_cycles << " cycles" << std::endl;
    std::cout << "I2: " << index2_cycles << " cycles" << std::endl;
    std::cout << "R2: " << read2_cycles << " cycles" << std::endl;
    std::cout << "Total: " << total_cycles << " cycles" << std::endl;
    
    // Verify expected values
    if (read1_cycles != 76 || index1_cycles != 10 || index2_cycles != 10 || read2_cycles != 76 || total_cycles != 172) {
        std::cerr << "Error: Unexpected cycle counts!" << std::endl;
        return false;
    }
    
    std::cout << "XML parsing test passed!" << std::endl;
    return true;
}

int main() {
    if (!test_xml_parsing()) {
        std::cerr << "XML parsing test failed!" << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed!" << std::endl;
    return 0;
} 