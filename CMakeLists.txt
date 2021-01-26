cmake_minimum_required(VERSION 3.12)

project(votca-tools LANGUAGES CXX)

set(PROJECT_VERSION "2022-dev")
string(REGEX REPLACE "[-.].*$" "" SOVERSION "${PROJECT_VERSION}")
if (NOT ${SOVERSION} MATCHES "^[0-9]+$")
  message(FATAL_ERROR "Could not determine SOVERSION (${SOVERSION}) from ${PROJECT_VERSION}")
endif (NOT ${SOVERSION} MATCHES "^[0-9]+$")
set(PROJECT_CONTACT "https://github.com/votca/tools/issues")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11


########################################################################
# User input options                                                   #
########################################################################
option(BUILD_SHARED_LIBS "Build shared libs" ON)
include(GNUInstallDirs)
include(FeatureSummary)

option(ENABLE_TESTING "Build and enable testing stuff" OFF)
add_feature_info(ENABLE_TESTING ENABLE_TESTING "Enable unit tests")
if(ENABLE_TESTING)
  enable_testing()
  include(FindPackageHandleStandardArgs) 
  find_program(VALGRIND_EXECUTABLE "valgrind")
  find_package_handle_standard_args(VALGRIND REQUIRED_VARS VALGRIND_EXECUTABLE)
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

find_package(Boost 1.53.0 REQUIRED COMPONENTS program_options filesystem system )
set_package_properties(Boost PROPERTIES TYPE REQUIRED PURPOSE "Extended C++ libraries")


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
set(MKL_THREAD_LAYER "Intel OpenMP" CACHE STRING "The thread layer to choose for MKL")
find_package(MKL)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set_package_properties(MKL PROPERTIES TYPE OPTIONAL PURPOSE "Enables accelerated performance with MKL")
endif (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)

find_package(Eigen3 3.3.0 NO_MODULE REQUIRED)
set_package_properties(Eigen3 PROPERTIES TYPE REQUIRED PURPOSE "C++ vector data structures")
message(STATUS "Found Eigen3: ${Eigen3_DIR}")

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

configure_file(${PROJECT_SOURCE_DIR}/CMakeModules/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
endif (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
