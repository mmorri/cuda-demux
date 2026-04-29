cmake_minimum_required(VERSION 3.16)

find_program(HOST_CXX_14 NAMES g++-14)
if(HOST_CXX_14 AND EXISTS "${HOST_CXX_14}" AND NOT DEFINED CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER "${HOST_CXX_14}" CACHE FILEPATH "" FORCE)
  set(CMAKE_CUDA_HOST_COMPILER "${HOST_CXX_14}")
  message(STATUS "Using host C++ compiler and CUDA host compiler: ${HOST_CXX_14}")
endif()

project(cuda-demux LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# Find required packages
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

# TinyXML2 detection: prefer system package, optional FetchContent
option(ENABLE_FETCH_TINYXML2 "Allow downloading TinyXML2 via FetchContent" OFF)

# Try config package names in common variants. Different distros export
# different target names (TinyXML2::TinyXML2 vs tinyxml2::tinyxml2 vs tinyxml2),
# so probe for the actual target after each find_package call.
set(TINYXML2_TARGET "")
find_package(TinyXML2 QUIET)
find_package(tinyxml2 QUIET)
foreach(_candidate TinyXML2::TinyXML2 tinyxml2::tinyxml2 tinyxml2)
  if(TARGET ${_candidate})
    set(TINYXML2_TARGET ${_candidate})
    break()
  endif()
endforeach()

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
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -allow-unsupported-compiler -Wno-deprecated-gpu-targets")

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
    CUDA_ARCHITECTURES "80;86;89"
)

# Link libraries
target_link_libraries(cuda-demux
    PRIVATE
    ZLIB::ZLIB
    ${TINYXML2_TARGET}
    OpenMP::OpenMP_CXX
    cudart
)

include(CTest)
if(BUILD_TESTING)
    add_executable(cuda-demux-unit-tests
        tests/test_demux.cpp
        src/fastq_writer.cpp
        src/common_adapters.cpp
    )
    target_include_directories(cuda-demux-unit-tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(cuda-demux-unit-tests PRIVATE ZLIB::ZLIB)
    add_test(NAME cuda-demux-unit-tests COMMAND cuda-demux-unit-tests)
endif()
