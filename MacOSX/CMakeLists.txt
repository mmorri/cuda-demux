cmake_minimum_required(VERSION 3.16)
project(cuda-demux-metal VERSION 1.0.0 LANGUAGES CXX OBJCXX)

# Check we're on macOS
if(NOT APPLE)
    message(FATAL_ERROR "This build configuration is only for macOS")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_OBJCXX_STANDARD 17)

# Find required frameworks
find_library(METAL_FRAMEWORK Metal)
find_library(METALPERFORMANCESHADERS_FRAMEWORK MetalPerformanceShaders)
find_library(FOUNDATION_FRAMEWORK Foundation)
find_library(COREIMAGE_FRAMEWORK CoreImage)

if(NOT METAL_FRAMEWORK)
    message(FATAL_ERROR "Metal framework not found")
endif()

# Find zlib
find_package(ZLIB REQUIRED)

# Find OpenMP (if available on macOS)
find_package(OpenMP)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${ZLIB_INCLUDE_DIRS}
)

# Compile Metal shaders to library
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/demux_kernels.metallib
    COMMAND xcrun -sdk macosx metal -c ${CMAKE_SOURCE_DIR}/shaders/demux_kernels.metal -o ${CMAKE_BINARY_DIR}/demux_kernels.air
    COMMAND xcrun -sdk macosx metallib ${CMAKE_BINARY_DIR}/demux_kernels.air -o ${CMAKE_BINARY_DIR}/demux_kernels.metallib
    DEPENDS ${CMAKE_SOURCE_DIR}/shaders/demux_kernels.metal
    COMMENT "Compiling Metal shaders"
)

add_custom_target(metal_shaders ALL DEPENDS ${CMAKE_BINARY_DIR}/demux_kernels.metallib)

# Source files
set(SOURCES
    src/bcl_parser_metal.mm
    src/demux_metal.mm
    src/main_metal.mm
    ../src/bcl_parser.cpp  # Reuse CPU implementation
    ../src/fastq_writer.cpp # Reuse FASTQ writer
    ../src/common_adapters.cpp # Reuse common code
)

# Create executable
add_executable(cuda-demux-metal ${SOURCES})

# Dependencies for shader compilation
add_dependencies(cuda-demux-metal metal_shaders)

# Link libraries
target_link_libraries(cuda-demux-metal
    ${METAL_FRAMEWORK}
    ${METALPERFORMANCESHADERS_FRAMEWORK}
    ${FOUNDATION_FRAMEWORK}
    ${ZLIB_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(cuda-demux-metal OpenMP::OpenMP_CXX)
endif()

# Copy shader library to output
add_custom_command(TARGET cuda-demux-metal POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_BINARY_DIR}/demux_kernels.metallib
    $<TARGET_FILE_DIR:cuda-demux-metal>/demux_kernels.metallib
)

# Installation
install(TARGETS cuda-demux-metal DESTINATION bin)
install(FILES ${CMAKE_BINARY_DIR}/demux_kernels.metallib DESTINATION lib/cuda-demux)

# Package configuration
set(CPACK_GENERATOR "DragNDrop")
set(CPACK_PACKAGE_NAME "CUDA-Demux-Metal")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_DESCRIPTION "GPU-accelerated Illumina demultiplexer for macOS using Metal")
set(CPACK_PACKAGE_VENDOR "NextGCon")
include(CPack)