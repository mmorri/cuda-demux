# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/mmorri/Desktop/test_cuda_demux/cuda-demux/build/_deps/tinyxml2-src")
  file(MAKE_DIRECTORY "/home/mmorri/Desktop/test_cuda_demux/cuda-demux/build/_deps/tinyxml2-src")
endif()
file(MAKE_DIRECTORY
  "/home/mmorri/Desktop/test_cuda_demux/cuda-demux/build/_deps/tinyxml2-build"
  "/home/mmorri/Desktop/test_cuda_demux/cuda-demux/build/_deps/tinyxml2-subbuild/tinyxml2-populate-prefix"
  "/home/mmorri/Desktop/test_cuda_demux/cuda-demux/build/_deps/tinyxml2-subbuild/tinyxml2-populate-prefix/tmp"
  "/home/mmorri/Desktop/test_cuda_demux/cuda-demux/build/_deps/tinyxml2-subbuild/tinyxml2-populate-prefix/src/tinyxml2-populate-stamp"
  "/home/mmorri/Desktop/test_cuda_demux/cuda-demux/build/_deps/tinyxml2-subbuild/tinyxml2-populate-prefix/src"
  "/home/mmorri/Desktop/test_cuda_demux/cuda-demux/build/_deps/tinyxml2-subbuild/tinyxml2-populate-prefix/src/tinyxml2-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/mmorri/Desktop/test_cuda_demux/cuda-demux/build/_deps/tinyxml2-subbuild/tinyxml2-populate-prefix/src/tinyxml2-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/mmorri/Desktop/test_cuda_demux/cuda-demux/build/_deps/tinyxml2-subbuild/tinyxml2-populate-prefix/src/tinyxml2-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
