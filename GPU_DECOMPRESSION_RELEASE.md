# CUDA-Demux GPU Decompression Release Notes

## Version 1.1.0-gpu

### Overview
This release introduces GPU-accelerated decompression for CBCL files, eliminating the CPU bottleneck that was causing 98% CPU usage during benchmarks.

### Problem Solved
The previous version (v1.0.1-cpu) showed high CPU usage (98%) because:
- CBCL files are compressed using GZIP/DEFLATE
- Decompression was happening sequentially on CPU using zlib
- Each tile and cycle required CPU decompression before GPU processing
- The GPU was idle waiting for decompressed data

### Solution Implemented
- **GPU Decompression Kernel**: Custom CUDA kernel for parallel GZIP/DEFLATE decompression
- **Batch Processing**: Multiple compressed blocks processed simultaneously on GPU
- **Zero CPU Overhead**: All decompression moved to GPU
- **Intelligent Fallback**: Automatic fallback to CPU if GPU decompression fails

### Performance Improvements
- **CPU Usage**: Reduced from 98% to minimal (only file I/O)
- **GPU Utilization**: Increased significantly
- **Throughput**: Faster processing of large CBCL files
- **Parallelization**: True parallel decompression of multiple tiles

### Technical Changes
1. **Removed Dependencies**:
   - OpenMP (CPU parallelization no longer needed)

2. **Added Components**:
   - `gpu_decompressor.cu`: CUDA kernel for GZIP decompression
   - `gpu_decompressor.h`: Header for GPU decompression functions

3. **Modified Files**:
   - `bcl_parser.cpp`: Updated to use GPU decompression
   - `CMakeLists.txt`: Removed OpenMP, added GPU decompressor
   - `main.cpp`: Added version information

### Usage

#### Default (GPU Decompression Enabled)
```bash
./cuda-demux --input /path/to/run --samplesheet sample.csv --output /path/to/output
```

#### Disable GPU Decompression (Use CPU Fallback)
```bash
export CUDA_DEMUX_NO_GPU_DECOMPRESS=1
./cuda-demux --input /path/to/run --samplesheet sample.csv --output /path/to/output
```

### Backward Compatibility
- The CPU version is preserved as tag `v1.0.1-cpu`
- GPU version includes CPU fallback for compatibility
- Environment variable control for decompression method

### Building from Source
```bash
git clone https://github.com/yourusername/cuda-demux.git
cd cuda-demux
git checkout v1.1.0-gpu
mkdir build && cd build
cmake .. -DENABLE_FETCH_TINYXML2=ON
make -j$(nproc)
```

### Requirements
- CUDA 11.0 or higher
- NVIDIA GPU with compute capability 8.0+ (Ampere or newer)
- Same requirements as v1.0.1 except OpenMP

### Known Limitations
- The GPU decompression kernel is optimized for GZIP format
- Very small files may not see significant speedup due to kernel launch overhead
- Fallback to CPU maintains compatibility but loses performance benefits

### Future Improvements
- Integration with nvCOMP for more compression formats
- Multi-GPU support for very large runs
- Further optimization of decompression kernel

### Benchmarking
To compare CPU vs GPU decompression:

```bash
# GPU decompression (default)
python3 scripts/benchmark_cuda_demux.py --base ~/Desktop/illumina --output ./benchmarks

# CPU decompression
export CUDA_DEMUX_NO_GPU_DECOMPRESS=1
python3 scripts/benchmark_cuda_demux.py --base ~/Desktop/illumina --output ./benchmarks-cpu
```

### Version Tags
- `v1.0.1-cpu`: Original version with CPU-based zlib decompression
- `v1.1.0-gpu`: New version with GPU-accelerated decompression