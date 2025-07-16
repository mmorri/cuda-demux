cmake_minimum_required(VERSION 3.16)
project(cuda-demux LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add the source files
add_executable(cuda-demux
    src/main.cpp
    src/demux.cu
    src/bcl_parser.cpp
    src/fastq_writer.cpp
)

# Set CUDA properties
set_target_properties(cuda-demux PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "52;60;61;70;75;80"
)
