# CUDA-Demux

**CUDA-Demux** is an open-source tool for demultiplexing Illumina BCL files into FASTQ format, optimized with CUDA for GPU acceleration.

## Features
- Fast barcode matching using CUDA.
- Multi-threaded BCL parsing.
- High-speed FASTQ generation.

## Repository
The code is hosted on GitHub at: [https://github.com/mmorri/cuda-demux](https://github.com/mmorri/cuda-demux).

## Installation
### Prerequisites
- CUDA Toolkit 11.0 or later
- CMake 3.16 or later
- A C++17-compatible compiler

### Build Instructions
1. Clone the repository:
   ```bash
   git clone https://github.com/mmorri/cuda-demux.git
   cd cuda-demux
