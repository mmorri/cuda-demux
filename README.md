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
- libdeflate (optional, `libdeflate-dev`): 2-3x faster gzip; zlib is used when absent

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
./cuda-demux --input <RUN_FOLDER> --samplesheet <SAMPLESHEET.CSV> --output <OUTPUT_FOLDER> [--gzip]

# Output / matching
  [--gzip]                  # write .fastq.gz (multi-member gzip, level 1 by default)
  [--gzip-level N]          # gzip level 1-9 (implies --gzip)
  [--barcode-mismatches N]  # mismatches allowed per index, 0-4 (default 1, like bcl-convert)

# Advanced GPU controls
  [--batch-size N]
  [--gpu-mem-fraction F]    # 0.05–0.95 fraction when cudaMemGetInfo unavailable
  [--device IDX]            # select CUDA device index
  [--verbose]               # per-file / per-batch diagnostics and phase timings
```

### Arguments

- `--input`: Path to the Illumina run folder containing the Data/Intensities/BaseCalls directory
- `--samplesheet`: Path to the CSV file with sample information and barcode mappings
- `--output`: Path to the directory where FASTQ files will be generated

### Output

One `<Sample_ID>_L<lane>_R1_001.fastq[.gz]` (and `_R2_` for paired runs) per sample and
lane, plus `undetermined_L<lane>_R{1,2}_001.fastq[.gz]`. Every sample in the sheet gets a
file even if no read matched it. Reads are written in cluster order and named
`@<Sample_ID>_<n>/1`, `@<Sample_ID>_<n>/2`. No-calls are emitted as `N` with quality `#`.
Quality scores use the bin table stored in each CBCL header.

### Index orientation

The i5 (Index2) orientation is taken from `RunInfo.xml` (`IsReverseComplement` on the
second index read) when present. Older runs without that attribute fall back to an
instrument-name heuristic from `RunParameters.xml`; `CUDA_DEMUX_I5_RC=0|1` overrides both.

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
MiSeq i100 Plus run, 32.2M passing-filter clusters, 2x151 + 2x10, 24 samples, `--gzip`,
RTX A4500 + 24-core/48-thread host: **8 s wall** (1.2 s CBCL ingest, 5 s demux + gzip),
14 GB RSS, 3.7 GB of output. The GPU decodes batch k+1 while the host formats and
compresses batch k on every core; each compressed chunk is an independent gzip member
(as produced by pigz / bcl-convert). gzip compression is the remaining bottleneck.

### Pipeline
1. `RunInfo.xml` gives the read structure and i5 orientation; lanes are processed one at
   a time (peak host memory is one lane: PF clusters x cycles bytes).
2. CBCL ingest (OpenMP over cycles): header bin table + `non_pf_excluded` flag, tiles
   mapped by id to the filter files, blocks inflated with libdeflate/zlib.
3. GPU: shared-memory transpose from cycle-major to read-major text, 2-bit barcode
   packing, popcount matching against the sample table in constant memory.
4. Host: reads grouped by sample, formatted and compressed in parallel, written in
   cluster order.

## Limitations

- The whole lane (PF clusters × cycles, one byte each) is held in host RAM while demuxing
- Requires sufficient GPU memory for barcode matching; batch size adapts to free memory
- Reads are named by ordinal, not by tile/x/y coordinates

## GPU Memory Sizing and Batching

The batch size (clusters per GPU pass) is derived from free device memory
(`cudaMemGetInfo`, or a 2 GiB assumption if unavailable) times a working fraction, capped
at 1M clusters; two batches are in flight at once. Overrides:

- `--batch-size N` / `CUDA_DEMUX_BATCH_SIZE`
- `--gpu-mem-fraction F` / `CUDA_DEMUX_MEM_FRACTION` (0.05–0.95, default 0.40)
- `--device N` / `CUDA_DEMUX_DEVICE`
- `--barcode-mismatches N` / `CUDA_DEMUX_MISMATCHES`
- `OMP_NUM_THREADS` for CBCL ingest and FASTQ compression threads

`--no-adaptive-probe` is accepted for compatibility but has no effect.

## Validation

`tests/reference_check.py` is an independent numpy CBCL decoder. It decodes a run straight
from the CBCL/filter files, applies the same matching rule, and compares a sample of reads
(every N-th cluster plus the first/last of the lane) against the tool's FASTQ output,
including per-sample read counts:

```bash
python3 tests/reference_check.py --run <RUN_FOLDER> --samplesheet <CSV> --output <OUT_DIR>
```

`tests/make_subset_run.py` cuts a few tiles out of a real run into a small, valid run
folder (default: the first tile of each surface), handy as a seconds-long regression
fixture for a given instrument layout:

```bash
python3 tests/make_subset_run.py --run <RUN_FOLDER> --out <SUBSET_FOLDER> --lanes 1
```

Unit tests (`cuda-demux-unit-tests`, built with `-DBUILD_TESTING=ON`) cover the FASTQ
writer, including multi-member gzip output and batch ordering.

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
