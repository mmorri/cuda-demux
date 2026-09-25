# CUDA-Demux Package Installation

This directory contains packaging configurations for various Linux distributions.

## Prerequisites

All packages require:
- NVIDIA GPU with compute capability 5.2 or higher
- CUDA Toolkit 11.0 or later
- zlib development libraries

## Package Types

### Arch Linux (AUR)

Build and install from the `arch/` directory:

```bash
cd arch/
makepkg -si
```

Or install from AUR (when published):
```bash
yay -S cuda-demux
```

### Debian/Ubuntu (.deb)

Install the .deb package:

```bash
sudo apt install ./cuda-demux_1.0.0-1_amd64.deb
```

Or build from source:
```bash
cd ..
dpkg-buildpackage -us -uc -b
sudo dpkg -i ../cuda-demux_*.deb
```

### Fedora/RHEL/CentOS (.rpm)

Install the RPM package:

```bash
sudo dnf install ./cuda-demux-1.0.0-1.x86_64.rpm
```

Or with yum:
```bash
sudo yum install ./cuda-demux-1.0.0-1.x86_64.rpm
```

## Building Packages

Use the provided build script to create packages:

```bash
# Build all packages
./build-packages.sh all

# Build specific package type
./build-packages.sh arch   # For Arch Linux
./build-packages.sh deb    # For Debian/Ubuntu
./build-packages.sh rpm    # For Fedora/RHEL
```

## Package Contents

All packages install:
- `/usr/bin/cuda-demux` - Main executable
- `/usr/share/doc/cuda-demux/README.md` - Documentation
- `/usr/share/licenses/cuda-demux/LICENSE` - License file

## Dependencies

### Runtime Dependencies
- CUDA runtime libraries (11.0+)
- zlib
- OpenMP runtime

### Build Dependencies
- CMake (>= 3.16)
- C++17 compatible compiler
- CUDA Toolkit
- zlib development headers
- OpenMP development files

## Troubleshooting

### CUDA Not Found
Ensure CUDA is properly installed and `nvcc` is in your PATH:
```bash
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
```

### Missing Dependencies
Install required development packages:

**Arch Linux:**
```bash
sudo pacman -S cuda zlib gcc openmp cmake
```

**Debian/Ubuntu:**
```bash
sudo apt install cuda-toolkit-11-0 zlib1g-dev libomp-dev cmake build-essential
```

**Fedora/RHEL:**
```bash
sudo dnf install cuda-toolkit zlib-devel libomp-devel cmake gcc-c++
```

## Testing Installation

After installation, test cuda-demux:

```bash
cuda-demux --help
```

## Uninstalling

**Arch Linux:**
```bash
sudo pacman -R cuda-demux
```

**Debian/Ubuntu:**
```bash
sudo apt remove cuda-demux
```

**Fedora/RHEL:**
```bash
sudo dnf remove cuda-demux
```