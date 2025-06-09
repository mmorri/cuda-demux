cmake_minimum_required(VERSION 3.16)
project(cuda-demux LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

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
)

# Set CUDA properties
set_target_properties(cuda-demux PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "52;60;61;70;75;80"
)

# Add test target
add_executable(test_demux
    tests/test_demux.cpp
    src/demux.cu
)
# Link to CUDA runtime if needed
find_package(CUDA REQUIRED)
target_link_libraries(cuda-demux PRIVATE CUDA::cudart)
target_link_libraries(test_demux PRIVATE CUDA::cudart)

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
