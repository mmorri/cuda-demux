cmake_minimum_required(VERSION 3.18)
project(cuda-demux LANGUAGES CXX CUDA)

enable_language(CUDA)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_ARCHITECTURES 75 80 86 CACHE STRING "CUDA architectures")

# Enable warnings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Add the source files
add_executable(cuda-demux
    src/main.cpp
    src/demux.cu
    src/demux.cpp
    src/bcl_parser.cpp
    src/fastq_writer.cpp
)
add_executable(test_demux
    tests/test_demux.cpp
    src/demux.cu
    src/demux.cpp
)
target_compile_features(cuda-demux PRIVATE cxx_std_17)
target_compile_features(test_demux PRIVATE cxx_std_17)

# Enable testing
include(CTest)
enable_testing()
add_test(NAME DemuxTest COMMAND test_demux)

# Doxygen documentation (optional)
find_package(Doxygen)
if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()
