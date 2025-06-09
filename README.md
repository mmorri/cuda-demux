# CUDA-Demux

**CUDA-Demux** is an open-source, GPU-accelerated tool for demultiplexing Illumina BCL files into FASTQ format. It uses CUDA for high-performance barcode matching and is designed for extensibility, robustness, and modern C++ best practices.

## Features
- Fast barcode matching using CUDA (GPU-accelerated)
- Modular, extensible C++17 codebase
- Robust error handling and logging
- High-speed FASTQ generation
- Unit and integration tests
- Doxygen documentation support

## Repository
The code is hosted on GitHub at: [https://github.com/mmorri/cuda-demux](https://github.com/mmorri/cuda-demux).

## Requirements
- CUDA Toolkit 11.0 or later
- CMake 3.16 or later
- A C++17-compatible compiler (GCC, Clang, or MSVC)
- NVIDIA GPU with compute capability 5.2 or higher

## Build Instructions
1. Clone the repository:
   ```bash
   git clone https://github.com/mmorri/cuda-demux.git
   cd cuda-demux
   ```
2. Create a build directory and navigate into it:
   ```bash
   mkdir build
   cd build
   ```
3. Run CMake to configure the project:
   ```bash
   cmake ..
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
./cuda-demux --input <BCL_FOLDER> --samplesheet <SAMPLESHEET.CSV> --output <OUTPUT_FOLDER>
```

### Arguments
- `--input`: Path to the directory containing .bcl files (raw sequencing data)
- `--samplesheet`: Path to the CSV file mapping barcodes to sample IDs
- `--output`: Path to the directory where FASTQ files will be generated

## Testing
To build and run tests:
```bash
make test_demux
./test_demux
```
Or run all tests with CTest:
```bash
ctest
```

## Documentation
If you have Doxygen installed, you can generate API documentation:
```bash
make doc
```
Documentation will be generated in the `doc/` directory inside your build folder.

## Troubleshooting
- Ensure your system has a compatible NVIDIA GPU and CUDA drivers installed
- If you encounter build errors, verify your CUDA and CMake versions
- For runtime errors, check GPU memory availability and input file integrity

## Contributing
Pull requests and issues are welcome! Please ensure code is well-documented and tested.

## License
This project is licensed under the MIT License.


