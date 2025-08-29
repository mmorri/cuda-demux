#include <iostream>
#include <tinyxml2.h>
#include <string>

using namespace tinyxml2;

int main() {
    // Test with the actual XML file
    XMLDocument doc;
    XMLError result = doc.LoadFile("test_runinfo.xml");
    
    if (result != XML_SUCCESS) {
        std::cerr << "XML parsing failed with error: " << doc.ErrorIDToName(result) << std::endl;
        return 1;
    }

    std::cout << "XML parsing successful!" << std::endl;

    // Test the same parsing logic as in bcl_parser.cpp
    XMLElement* run_info_element = doc.FirstChildElement("RunInfo");
    if (!run_info_element) {
        std::cerr << "Error: <RunInfo> element not found in RunInfo.xml." << std::endl;
        return 1;
    }

    XMLElement* run_element = run_info_element->FirstChildElement("Run");
    if (!run_element) {
        std::cerr << "Error: <Run> element not found in RunInfo.xml." << std::endl;
        return 1;
    }

    XMLElement* reads_element = run_element->FirstChildElement("Reads");
    if (!reads_element) {
        std::cerr << "Error: <Reads> element not found in RunInfo.xml." << std::endl;
        return 1;
    }

    std::cout << "Found RunInfo, Run, and Reads elements successfully!" << std::endl;

    // Parse reads like in the original code
    int read1_cycles = 0;
    int index1_cycles = 0;
    int index2_cycles = 0;
    int read2_cycles = 0;
    int total_cycles = 0;
    int read_count = 0;

    std::cout << "Parsing Read elements..." << std::endl;
    for (XMLElement* read_elem = reads_element->FirstChildElement("Read"); 
         read_elem != nullptr; 
         read_elem = read_elem->NextSiblingElement("Read")) {
        
        read_count++;
        
        const char* num_cycles_attr = read_elem->Attribute("NumCycles");
        const char* is_indexed_attr = read_elem->Attribute("IsIndexedRead");
        
        if (!num_cycles_attr || !is_indexed_attr) {
            std::cerr << "Error: Missing required attributes in Read element " << read_count << std::endl;
            return 1;
        }
        
        int num_cycles = std::stoi(num_cycles_attr);
        bool is_indexed = (std::string(is_indexed_attr) == "Y");
        
        std::cout << "Read " << read_count << ": NumCycles=" << num_cycles << ", IsIndexed=" << (is_indexed ? "Y" : "N") << std::endl;
        
        if (!is_indexed) {
            if (read1_cycles == 0) { 
                read1_cycles = num_cycles; 
                std::cout << "  -> Assigned to Read1" << std::endl;
            }
            else { 
                read2_cycles = num_cycles; 
                std::cout << "  -> Assigned to Read2" << std::endl;
            }
        } else {
            if (index1_cycles == 0) { 
                index1_cycles = num_cycles; 
                std::cout << "  -> Assigned to Index1" << std::endl;
            }
            else { 
                index2_cycles = num_cycles; 
                std::cout << "  -> Assigned to Index2" << std::endl;
            }
        }
        total_cycles += num_cycles;
    }

    std::cout << "Parsed " << read_count << " Read elements" << std::endl;
    std::cout << "\nFinal Run Structure:" << std::endl;
    std::cout << "R1: " << read1_cycles << " cycles" << std::endl;
    std::cout << "I1: " << index1_cycles << " cycles" << std::endl;
    std::cout << "I2: " << index2_cycles << " cycles" << std::endl;
    std::cout << "R2: " << read2_cycles << " cycles" << std::endl;
    std::cout << "Total: " << total_cycles << " cycles" << std::endl;

    // Verify expected values
    if (read1_cycles != 76 || index1_cycles != 10 || index2_cycles != 10 || read2_cycles != 76 || total_cycles != 172) {
        std::cerr << "Error: Unexpected cycle counts!" << std::endl;
        return 1;
    }

    std::cout << "XML parsing test passed!" << std::endl;
    return 0;
} 