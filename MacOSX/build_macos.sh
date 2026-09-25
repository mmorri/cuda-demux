#!/bin/bash

# Build script for macOS Metal version of cuda-demux

set -e

echo "================================"
echo "CUDA-Demux Metal Build for macOS"
echo "================================"

# Check if we're on macOS
if [[ "$OSTYPE" != "darwin"* ]]; then
    echo "Error: This script is for macOS only"
    exit 1
fi

# Check for Metal support
if ! xcrun --sdk macosx metal --version &>/dev/null; then
    echo "Error: Metal compiler not found. Please install Xcode Command Line Tools"
    echo "Run: xcode-select --install"
    exit 1
fi

# Parse arguments
BUILD_TYPE="Release"
INSTALL_PREFIX="/usr/local"

while [[ $# -gt 0 ]]; do
    case $1 in
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        --prefix)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --debug          Build in debug mode"
            echo "  --prefix <path>  Installation prefix (default: /usr/local)"
            echo "  --help           Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create build directory
BUILD_DIR="build-metal"
mkdir -p $BUILD_DIR
cd $BUILD_DIR

echo ""
echo "Configuration:"
echo "  Build Type: $BUILD_TYPE"
echo "  Install Prefix: $INSTALL_PREFIX"
echo ""

# Configure with CMake
echo "Configuring with CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
    -DCMAKE_OSX_DEPLOYMENT_TARGET=11.0

# Build
echo ""
echo "Building..."
cmake --build . -j$(sysctl -n hw.ncpu)

# Run tests if they exist
if [[ -f "test/test_metal" ]]; then
    echo ""
    echo "Running tests..."
    ./test/test_metal
fi

echo ""
echo "Build complete!"
echo ""
echo "To install, run:"
echo "  cd $BUILD_DIR && sudo cmake --install ."
echo ""
echo "To create a macOS app bundle:"
echo "  cd $BUILD_DIR && cpack -G DragNDrop"