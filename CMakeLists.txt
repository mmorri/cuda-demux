cmake_minimum_required(VERSION 3.16)
project(cuda-demux LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CUDA_STANDARD 17)

# Find required packages
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)
find_package(tinyxml2 QUIET)

# If tinyxml2 is not found, download and build it
if(NOT tinyxml2_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    tinyxml2
    GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
    GIT_TAG 9.0.0
  )
  FetchContent_MakeAvailable(tinyxml2)
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add CUDA BCL parser file
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")

# Add the source files
add_executable(cuda-demux
    src/main.cpp
    src/demux.cu
    src/bcl_parser.cpp
    src/bcl_parser_cuda.cu
    src/fastq_writer.cpp
    src/common_adapters.cpp
)

# Set CUDA properties
set_target_properties(cuda-demux PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "52;60;61;70;75;80"
)

# Link libraries
target_link_libraries(cuda-demux
    PRIVATE
    ZLIB::ZLIB
    tinyxml2
    OpenMP::OpenMP_CXX
    cudart
)
