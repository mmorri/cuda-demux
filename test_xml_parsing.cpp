#include <iostream>
#include <tinyxml2.h>
#include <string>

using namespace tinyxml2;

int main() {
    // Test with the provided XML content
    const char* xml_content = R"(
<?xml version="1.0"?>
<RunInfo Version="6">
	<Run Id="240227_VH00103_437_2222FFYNX" Number="437">
		<Flowcell>2222FFYNX</Flowcell>
		<Instrument>VH00103</Instrument>
		<Date>2024-02-27T20:40:35Z</Date>
		<Reads>
			<Read Number="1" NumCycles="76" IsIndexedRead="N" IsReverseComplement="N"/>
			<Read Number="2" NumCycles="10" IsIndexedRead="Y" IsReverseComplement="N"/>
			<Read Number="3" NumCycles="10" IsIndexedRead="Y" IsReverseComplement="Y"/>
			<Read Number="4" NumCycles="76" IsIndexedRead="N" IsReverseComplement="N"/>
		</Reads>
		<FlowcellLayout LaneCount="2" SurfaceCount="2" SwathCount="6" TileCount="16">
			<TileSet TileNamingConvention="FourDigit">
				<Tiles>
					<Tile>1_1101</Tile>
					<Tile>1_1102</Tile>
				</Tiles>
			</TileSet>
		</FlowcellLayout>
		<ImageDimensions Width="8208" Height="5541"/>
		<ImageChannels>
			<Name>green</Name>
			<Name>blue</Name>
		</ImageChannels>
	</Run>
</RunInfo>
)";

    XMLDocument doc;
    XMLError result = doc.Parse(xml_content);
    
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

    for (XMLElement* read_elem = reads_element->FirstChildElement("Read"); 
         read_elem != nullptr; 
         read_elem = read_elem->NextSiblingElement("Read")) {
        
        int num_cycles = std::stoi(read_elem->Attribute("NumCycles"));
        bool is_indexed = (std::string(read_elem->Attribute("IsIndexedRead")) == "Y");
        
        std::cout << "Read: NumCycles=" << num_cycles << ", IsIndexed=" << (is_indexed ? "Y" : "N") << std::endl;
        
        if (!is_indexed) {
            if (read1_cycles == 0) { 
                read1_cycles = num_cycles; 
                std::cout << "  -> Read1 cycles: " << num_cycles << std::endl;
            }
            else { 
                read2_cycles = num_cycles; 
                std::cout << "  -> Read2 cycles: " << num_cycles << std::endl;
            }
        } else {
            if (index1_cycles == 0) { 
                index1_cycles = num_cycles; 
                std::cout << "  -> Index1 cycles: " << num_cycles << std::endl;
            }
            else { 
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

    return 0;
} 