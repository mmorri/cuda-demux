# CUDA-Demux

**CUDA-Demux** is a high-performance, GPU-accelerated tool for demultiplexing Illumina sequencing data from BCL/CBCL format to FASTQ files. It leverages CUDA for fast barcode matching and parallel processing.

## Features
- **GPU Acceleration**: Fast barcode matching using CUDA kernels
- **CBCL Support**: Native support for compressed BCL (CBCL) files from NovaSeq and NextSeq platforms
- **Multi-threaded**: Parallel BCL/CBCL file parsing and decompression
- **Dual Indexing**: Support for single and dual-indexed libraries
- **High Performance**: Process billions of reads efficiently

## Repository
The code is hosted on GitHub at: [https://github.com/mmorri/cuda-demux](https://github.com/mmorri/cuda-demux).

## Installation
### Prerequisites
- CUDA Toolkit 11.0 or later
- CMake 3.16 or later
- A C++17-compatible compiler (GCC 7+ or Clang 5+)
- NVIDIA GPU with compute capability 5.2 or higher
- zlib development libraries
- OpenMP support

### Build Instructions
1. Clone the repository:
   ```bash
   git clone https://github.com/mmorri/cuda-demux.git
   cd cuda-demux
   ```
2. Create a build directory:
   ```bash
   mkdir build
   cd build
   ```
3. Configure the project with CMake:
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Release
   ```
4. Compile the tool:
   ```bash
   make -j$(nproc)
   ```
5. Verify the binary is created:
   ```bash
   ls cuda-demux
   ```

## Usage

```bash
./cuda-demux --input <RUN_FOLDER> --samplesheet <SAMPLESHEET.CSV> --output <OUTPUT_FOLDER>
```

### Arguments

- `--input`: Path to the Illumina run folder containing the Data/Intensities/BaseCalls directory
- `--samplesheet`: Path to the CSV file with sample information and barcode mappings
- `--output`: Path to the directory where FASTQ files will be generated

### Example

```bash
./cuda-demux \
  --input /path/to/NovaSeqX_Run \
  --samplesheet /path/to/SampleSheet.csv \
  --output /path/to/output
```

## Input Requirements

### Run Folder Structure
The tool expects a standard Illumina run folder structure:
```
Run_Folder/
├── RunInfo.xml
├── RunParameters.xml
├── SampleSheet.csv
└── Data/
    └── Intensities/
        └── BaseCalls/
            └── L001/
                ├── C1.1/
                │   ├── L001_1.cbcl
                │   └── L001_2.cbcl
                ├── C2.1/
                └── ...
```

### Sample Sheet Format
The sample sheet should follow the Illumina format with sections for:
- `[Header]` - Run metadata
- `[Reads]` - Read structure
- `[BCLConvert_Data]` - Sample barcode mappings

Example:
```csv
[Header]
FileFormatVersion,2
RunName,MyRun
InstrumentPlatform,NextSeq1k2k

[Reads]
Read1Cycles,151
Read2Cycles,151
Index1Cycles,10
Index2Cycles,10

[BCLConvert_Data]
Sample_ID,Index,Index2
Sample1,ATCGATCGAT,TAGCTAGCTA
Sample2,GCTAGCTAGC,CGATCGATCG
```

## Performance

The tool is optimized for:
- Large-scale sequencing runs (tested with >1 billion clusters)
- NovaSeq, NextSeq, and other Illumina platforms using CBCL format
- Multi-GPU systems (future enhancement)

### Benchmarks
- Processes ~10 million reads in minutes on a modern GPU
- Scales linearly with number of reads
- Memory usage depends on batch size and number of samples

## Limitations

- Currently processes reads in batches (default: 10M reads for testing)
- Requires sufficient GPU memory for barcode matching
- Output files are uncompressed FASTQ (gzip compression planned)

## Troubleshooting

1. **CUDA errors**: Ensure your GPU driver and CUDA toolkit are properly installed
2. **Memory issues**: Reduce batch size or ensure sufficient GPU memory
3. **File not found**: Verify the run folder structure matches Illumina standards

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Uses tinyxml2 for XML parsing
- CUDA toolkit for GPU acceleration
- zlib for CBCL decompression





