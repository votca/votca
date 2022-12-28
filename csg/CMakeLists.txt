cmake_minimum_required(VERSION 3.13)

project(votca-csg LANGUAGES CXX)

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
  set_package_properties(Git PROPERTIES TYPE OPTIONAL PURPOSE "Generated version for development version")
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)

# iie methods needs python 3.5
find_package(Python 3.5 REQUIRED COMPONENTS Interpreter)
set_package_properties(Python PROPERTIES TYPE OPTIONAL PURPOSE "Used for csg scripts")
option(BUILD_MANPAGES "Build manpages" ${Python_FOUND})
add_feature_info(BUILD_MANPAGES BUILD_MANPAGES "Build manpages (disable for cross-compile)")

find_package(Boost 1.71.0 REQUIRED COMPONENTS program_options filesystem system regex)
find_package(Eigen3 3.3.0 NO_MODULE REQUIRED)
find_package(Threads REQUIRED)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set_package_properties(Threads PROPERTIES TYPE REQUIRED PURPOSE "Used for thread parallelization")
  set_package_properties(Boost PROPERTIES TYPE REQUIRED PURPOSE "Extended C++ libraries")
  set_package_properties(Eigen3 PROPERTIES TYPE REQUIRED PURPOSE "C++ vector data structures")
endif()
message(STATUS "Found Eigen3: ${Eigen3_DIR}")

find_program(GMX_EXECUTABLE NAMES gmx_d gmx)
find_package_handle_standard_args(GMX REQUIRED_VARS GMX_EXECUTABLE)

if(ENABLE_TESTING)
  set(INTEGRATIONTEST_TOLERANCE "5e-5" CACHE STRING "Tolerance for the regression tests")
  mark_as_advanced(INTEGRATIONTEST_TOLERANCE)
  find_package(UnixCommands)
  set_package_properties(UnixCommands PROPERTIES TYPE OPTIONAL PURPOSE "Used in testing")
endif(ENABLE_TESTING)

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFileCXX)
foreach(HEADER algorithm cassert fstream iomanip iostream list map numeric sstream stdexcept string vector cstdlib)
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
if(NOT TARGET manpages)
  add_custom_target(manpages ALL)
endif()
if(VOTCA_SPHINX_DIR)
  add_custom_target(doc-csg)
  file(MAKE_DIRECTORY ${VOTCA_SPHINX_DIR}/csg)
endif()
add_subdirectory(src)
add_subdirectory(scripts)
add_subdirectory(include/votca/csg)
add_subdirectory(share)

# needs to happen after src and scripts subdir
get_property(VOTCA_BINARIES TARGET votca_csg PROPERTY BINARIES)
list(SORT VOTCA_BINARIES)
configure_file(src/libcsg/VOTCA_CSGConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/src/libcsg/VOTCA_CSGConfig.cmake" @ONLY)
list(JOIN VOTCA_BINARIES " " VOTCA_BINARIES)
configure_file(scripts/csg-completion.bash.in "${CMAKE_CURRENT_BINARY_DIR}/scripts/csg-completion.bash" @ONLY)
