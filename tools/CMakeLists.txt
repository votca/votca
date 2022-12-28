cmake_minimum_required(VERSION 3.13)

project(votca-tools LANGUAGES CXX)

set(PROJECT_VERSION "2023-dev")

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(CheckCXXCompilerFlag)

########################################################################
# User input options                                                   #
########################################################################
option(BUILD_SHARED_LIBS "Build shared libs" ON)
include(GNUInstallDirs)

option(ENABLE_TESTING "Build and enable testing stuff" OFF)
add_feature_info(ENABLE_TESTING ENABLE_TESTING "Enable unit tests")
if(ENABLE_TESTING)
  enable_testing()
endif(ENABLE_TESTING)

########################################################################
#Find external packages
########################################################################
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
  set_package_properties(Git PROPERTIES TYPE OPTIONAL PURPOSE "Generates version for development version")
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)

find_package(Threads REQUIRED)
set_package_properties(Threads PROPERTIES TYPE REQUIRED PURPOSE "Used for thread parallelization")

find_package(Python 3 COMPONENTS Interpreter)
set_package_properties(Python PROPERTIES TYPE OPTIONAL PURPOSE "Used to build documentation files")
option(BUILD_MANPAGES "Build manpages (might lead to problem on system without rpath" ${Python_FOUND})
add_feature_info(BUILD_MANPAGES BUILD_MANPAGES "Build manpages (disable for cross-compile)")
#define this target here, so that individual man pages can append to it.
add_custom_target(manpages ALL)
if(VOTCA_SPHINX_DIR)
  file(MAKE_DIRECTORY ${VOTCA_SPHINX_DIR}/tools)
  add_custom_target(doc-tools)
endif()

########################################################################
# Checks what linear algebra packages are installed                    #
########################################################################
find_package(MKL CONFIG)
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

set(MATH_LIBRARIES "m" CACHE STRING "math library")
mark_as_advanced( MATH_LIBRARIES )

######################################
# Include the following subdirectory #
######################################
add_subdirectory(scripts)
add_subdirectory(src)
add_subdirectory(include/votca/tools)
add_subdirectory(share/man)
add_subdirectory(share/doc)
