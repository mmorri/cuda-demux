# AGENTS.md

## Cursor Cloud specific instructions

### Product
`cuda-demux` is a **CLI binary** (C++17 + CUDA), not a web service. There is no local app server, database, or login flow. End-to-end work means: configure with CMake, build the binary, then run it against an Illumina run folder.

### Build (standard commands)
See `README.md` and `.github/workflows/build-test.yml`. Preferred configure flags in this environment:

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=gcc-12 \
  -DCMAKE_CXX_COMPILER=g++-12 \
  -DCMAKE_CUDA_HOST_COMPILER=g++-12
make -j$(nproc)
```

System tinyxml2 is installed (`libtinyxml2-dev`). If it is missing, use `-DENABLE_FETCH_TINYXML2=ON` (as CI does).

### Compiler gotchas
- Ubuntu’s `c++` alternative may point at clang; CUDA builds should use **GCC**. Prefer `g++-12` as the CUDA host compiler with the apt `nvidia-cuda-toolkit` (12.0) package — newer GCC hosts often fail with older nvcc.
- Shell profile in this VM sets `CC=gcc-12`, `CXX=g++-12`, `CUDAHOSTCXX=g++-12`.

### Lint / test
- No ESLint/clang-tidy target is wired in CMake. Treat successful `cmake` + `make` as the compile check.
- `tests/test_demux.cpp` is a stub and is **not** a CMake target. Compile manually if needed:
  `g++-12 -std=c++17 -I include tests/test_demux.cpp -o build/test_demux`
- CI’s “basic test” is `build/cuda-demux --help` / `--version`.

### Run / GPU constraint
- This Cloud Agent VM has **no NVIDIA GPU**. `demux()` and BCL CUDA decode require a CUDA device; without one the process exits early with “No CUDA-capable GPU found”.
- For real demux runs, use a remote GPU host over SSH (or a small HTTP API on that host). Expected secrets when configured:
  - `GPU_SSH_HOST` (e.g. `98.210.223.214`)
  - `GPU_SSH_USER`
  - `GPU_SSH_PASSWORD` (or prefer key-based auth via `GPU_SSH_PRIVATE_KEY`)
- Do **not** commit passwords or private keys into the repo.

### Smoke commands (no GPU)
```bash
./build/cuda-demux --version
./build/cuda-demux   # prints usage; exit 1 is expected
```
