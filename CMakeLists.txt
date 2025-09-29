cmake_minimum_required(VERSION 3.16)

find_program(HOST_CXX_14 NAMES g++-14)
if(HOST_CXX_14)
  set(CMAKE_CXX_COMPILER "${HOST_CXX_14}" CACHE FILEPATH "" FORCE)
  set(CMAKE_CUDA_HOST_COMPILER "${HOST_CXX_14}")
  message(STATUS "Using host C++ compiler and CUDA host compiler: ${HOST_CXX_14}")
endif()

# Allow unsupported compiler for CUDA and set C++ standard
set(CMAKE_CUDA_FLAGS "-allow-unsupported-compiler")

project(cuda-demux LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# Find required packages
find_package(ZLIB REQUIRED)
find_package(CUDAToolkit REQUIRED)

# TinyXML2 detection: prefer system package, optional FetchContent
option(ENABLE_FETCH_TINYXML2 "Allow downloading TinyXML2 via FetchContent" OFF)

# Try config package names in common variants
set(TINYXML2_TARGET "")
find_package(TinyXML2 QUIET)
if(TinyXML2_FOUND)
  set(TINYXML2_TARGET TinyXML2::TinyXML2)
else()
  find_package(tinyxml2 QUIET)
  if(tinyxml2_FOUND)
    # Some builds export target 'tinyxml2'
    if(TARGET tinyxml2::tinyxml2)
      set(TINYXML2_TARGET tinyxml2::tinyxml2)
    elseif(TARGET tinyxml2)
      set(TINYXML2_TARGET tinyxml2)
    endif()
  endif()
endif()

if(NOT TINYXML2_TARGET)
  if(ENABLE_FETCH_TINYXML2)
    include(FetchContent)
    FetchContent_Declare(
      tinyxml2
      GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
      GIT_TAG 9.0.0
    )
    FetchContent_MakeAvailable(tinyxml2)
    set(TINYXML2_TARGET tinyxml2)
  else()
    message(FATAL_ERROR "tinyxml2 not found. Install libtinyxml2-dev or re-run with -DENABLE_FETCH_TINYXML2=ON to fetch.")
  endif()
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# CUDA flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr -allow-unsupported-compiler -Wno-deprecated-gpu-targets")

# Add the source files
add_executable(cuda-demux
    src/main.cpp
    src/demux.cu
    src/bcl_parser.cpp
    src/bcl_parser_cuda.cu
    src/gpu_decompressor.cu
    src/fastq_writer.cpp
    src/common_adapters.cpp
    src/memory_utils.cpp
)

# Set CUDA properties
set_target_properties(cuda-demux PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "80;86;89"
)

# Link libraries
target_link_libraries(cuda-demux
    PRIVATE
    ZLIB::ZLIB
    tinyxml2
    cudart
)
