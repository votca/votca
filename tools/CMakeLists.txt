cmake_minimum_required(VERSION 3.13)

project(votca-tools LANGUAGES CXX)

set(PROJECT_VERSION "2023-dev")

########################################################################
#Find external packages
########################################################################
if(VOTCA_SPHINX_DIR)
  file(MAKE_DIRECTORY ${VOTCA_SPHINX_DIR}/tools)
  add_custom_target(doc-tools)
endif()

########################################################################
# Checks what linear algebra packages are installed                    #
########################################################################

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(MKL_THREAD_LAYER "Intel OpenMP" CACHE STRING "The thread layer to choose for MKL")
else()
  set(MKL_THREAD_LAYER "GNU OpenMP" CACHE STRING "The thread layer to choose for MKL")
endif()
find_package(MKL)
set_package_properties(MKL PROPERTIES TYPE OPTIONAL PURPOSE "Enables accelerated performance with MKL")
# Replace the below with CMAKE_REQUIRE_FIND_PACKAGE_MKL, when we switch to
# CMake-3.22+
option(REQUIRE_MKL "Require MKL" OFF)
if(REQUIRE_MKL)
  find_package(MKL REQUIRED)
endif()

#user defined reductions are buggy for <= clang-9 (see https://bugs.llvm.org/show_bug.cgi?id=44134)
if (CMAKE_CXX_COMPILER MATCHES "clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
  # in case OpenMP was detected somewhere else
  unset(OpenMP_FOUND)
else()
  find_package(OpenMP 4.5)
  set_package_properties(OpenMP PROPERTIES TYPE RECOMMENDED PURPOSE "Used for thread parallelization in tools")
endif()

if(MKL_FOUND)
  # MKL has its internal fftw3 implementation
  unset(FFTW3_FOUND)
else()
  find_package(FFTW3)
  set_package_properties(FFTW3 PROPERTIES TYPE RECOMMENDED PURPOSE "Enables FFT functionality")
endif()

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFileCXX)
foreach(HEADER cassert cmath fstream functional iostream limits list map ostream sstream stack stdexcept string vector)
  check_include_file_cxx(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

######################################
# Include the following subdirectory #
######################################
add_subdirectory(scripts)
add_subdirectory(src)
add_subdirectory(include/votca/tools)
add_subdirectory(share/man)
add_subdirectory(share/doc)
