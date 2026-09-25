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

# libdeflate (optional): 2-3x faster than zlib for both the gzip output and
# CBCL inflate, and gzip is the bottleneck once the GPU path is fast. Static
# linking is preferred so the binary gains no extra runtime dependency, which
# matters for the released .deb/.rpm: their dependency lists name only libc,
# libgcc, libstdc++, zlib and libgomp.
option(ENABLE_FETCH_LIBDEFLATE "Download and statically link libdeflate" OFF)

set(HAVE_LIBDEFLATE OFF)
set(LIBDEFLATE_TARGET "")

if(ENABLE_FETCH_LIBDEFLATE)
  include(FetchContent)
  set(LIBDEFLATE_BUILD_SHARED_LIB OFF CACHE BOOL "" FORCE)
  set(LIBDEFLATE_BUILD_GZIP OFF CACHE BOOL "" FORCE)
  set(LIBDEFLATE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(
    libdeflate
    GIT_REPOSITORY https://github.com/ebiggers/libdeflate.git
    GIT_TAG v1.23
  )
  FetchContent_MakeAvailable(libdeflate)
  set(LIBDEFLATE_TARGET libdeflate::libdeflate_static)
  set(HAVE_LIBDEFLATE ON)
  message(STATUS "Using libdeflate: fetched v1.23, linked statically")
else()
  # Prefer the static archive over the shared object for the same reason.
  find_path(LIBDEFLATE_INCLUDE_DIR libdeflate.h)
  find_library(LIBDEFLATE_LIBRARY NAMES libdeflate.a deflate)
  if(LIBDEFLATE_INCLUDE_DIR AND LIBDEFLATE_LIBRARY)
    add_library(cuda_demux_libdeflate INTERFACE)
    target_include_directories(cuda_demux_libdeflate INTERFACE ${LIBDEFLATE_INCLUDE_DIR})
    target_link_libraries(cuda_demux_libdeflate INTERFACE ${LIBDEFLATE_LIBRARY})
    set(LIBDEFLATE_TARGET cuda_demux_libdeflate)
    set(HAVE_LIBDEFLATE ON)
    message(STATUS "Using libdeflate: ${LIBDEFLATE_LIBRARY}")
  else()
    message(STATUS "libdeflate not found; gzip output will use zlib. Install "
                   "libdeflate-dev, or configure with -DENABLE_FETCH_LIBDEFLATE=ON.")
  endif()
endif()

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
# nvcc does not forward OpenMP flags to its host compiler by itself; the .cu
# files use OpenMP for host-side staging.
target_compile_options(cuda-demux PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>
)

# Link libraries
target_link_libraries(cuda-demux
    PRIVATE
    ZLIB::ZLIB
    ${TINYXML2_TARGET}
    OpenMP::OpenMP_CXX
    cudart
)
if(HAVE_LIBDEFLATE)
  target_compile_definitions(cuda-demux PRIVATE CUDA_DEMUX_HAVE_LIBDEFLATE=1)
  target_link_libraries(cuda-demux PRIVATE ${LIBDEFLATE_TARGET})
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(cuda-demux-unit-tests
        tests/test_demux.cpp
        src/fastq_writer.cpp
        src/common_adapters.cpp
    )
    target_include_directories(cuda-demux-unit-tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
    # The tests check with a macro that is not stripped by -DNDEBUG in Release.
    target_compile_definitions(cuda-demux-unit-tests PRIVATE CUDA_DEMUX_TESTS=1)
    target_link_libraries(cuda-demux-unit-tests PRIVATE ZLIB::ZLIB OpenMP::OpenMP_CXX)
    if(HAVE_LIBDEFLATE)
      target_compile_definitions(cuda-demux-unit-tests PRIVATE CUDA_DEMUX_HAVE_LIBDEFLATE=1)
      target_link_libraries(cuda-demux-unit-tests PRIVATE ${LIBDEFLATE_TARGET})
    endif()
    add_test(NAME cuda-demux-unit-tests COMMAND cuda-demux-unit-tests)
endif()
