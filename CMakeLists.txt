cmake_minimum_required(VERSION 3.10)
project(votca-xtp LANGUAGES CXX)

set(PROJECT_VERSION "2021-dev")
string(REGEX REPLACE "[-.].*$" "" SOVERSION "${PROJECT_VERSION}")
if (NOT ${SOVERSION} MATCHES "^[0-9]+$")
  message(FATAL_ERROR "Could not determine SOVERSION (${SOVERSION}) from ${PROJECT_VERSION}")
endif (NOT ${SOVERSION} MATCHES "^[0-9]+$")
set(PROJECT_CONTACT "https://github.com/votca/xtp/issues")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_definitions(-DDEBUG)
endif()

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(FeatureSummary)
include(CheckCXXCompilerFlag)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

#openmp with clang smaller version 10 and user defined reductions is buggy https://bugs.llvm.org/show_bug.cgi?id=44134
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)
  find_package(OpenMP 4.5)
  set_package_properties(OpenMP PROPERTIES TYPE RECOMMENDED PURPOSE "Used for thread parallelization in xtp")
else()
  # in case OpenMP was detected somewhere else
  unset(OpenMP_FOUND)
endif()

########################################################################
# User input options                                                   #
########################################################################

find_package(PkgConfig REQUIRED)
set_package_properties(PkgConfig PROPERTIES TYPE RECOMMENDED PURPOSE "Used to detect libxc and libint package")
pkg_check_modules(LIBXC REQUIRED IMPORTED_TARGET libxc)
pkg_check_modules(LIBINT REQUIRED IMPORTED_TARGET libint2)

find_package(HDF5 1.8 REQUIRED COMPONENTS CXX)
set_package_properties(HDF5 PROPERTIES TYPE REQUIRED PURPOSE "Used to read/write HDF5 data files")

if(HDF5_VERSION VERSION_GREATER 1.8)
  message(WARNING "HDF5 will be used such that it is compatible with version 1.8.")
endif()

# https://github.com/votca/xtp/issues/436, hdf5-1.10.4 generates a implicitly-declared operator warning
if(HDF5_VERSION VERSION_GREATER_EQUAL 1.10.4 AND HDF5_VERSION VERSION_LESS_EQUAL 1.10.6)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-Wno-deprecated-copy" COMPILER_SUPPORTS_WNO_DEPRECATED_COPY)
  if(COMPILER_SUPPORTS_WNO_DEPRECATED_COPY)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-copy")
  endif()
endif()

option(BUILD_SHARED_LIBS "Build shared libs" ON)
include(GNUInstallDirs)

option(USE_CUDA "Use a Nvidia GPU to speed up computations" OFF)
add_feature_info(USE_CUDA USE_CUDA "Use CUDA for GPU acceleration")

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

find_package(Boost 1.57.0 REQUIRED COMPONENTS program_options filesystem system timer)
find_package(Eigen3 3.3.0 REQUIRED NO_MODULE)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set_package_properties(Boost PROPERTIES TYPE REQUIRED PURPOSE "Extended C++ libraries")
  set_package_properties(Eigen3 PROPERTIES TYPE REQUIRED PURPOSE "C++ vector data structures")
endif()
message(STATUS "Found Eigen3: ${Eigen3_DIR}")

option(BUILD_MANPAGES "Build manpages (might lead to problem on system without rpath" OFF)
add_feature_info(BUILD_MANPAGES BUILD_MANPAGES "Build manpages (disable for cross-compile)")
#define this target here, so that individual man pages can append to it.
if(NOT TARGET manpages)
  add_custom_target(manpages ALL)
endif()

find_package(VOTCA_TOOLS REQUIRED NO_MODULE)
set_package_properties(VOTCA_TOOLS PROPERTIES TYPE REQUIRED PURPOSE "Votca base library")
if(DEFINED VOTCA_TOOLS_VERSION AND NOT VOTCA_TOOLS_VERSION STREQUAL PROJECT_VERSION)
  message(FATAL_ERROR "Incompatible VOTCA Tools version found (needed ${PROJECT_VERSION}, found ${VOTCA_TOOLS_VERSION})")
endif()
message(STATUS "Found VOTCA_TOOLS: ${VOTCA_TOOLS_DIR}")
find_package(VOTCA_CSG REQUIRED NO_MODULE)
set_package_properties(VOTCA_TOOLS PROPERTIES TYPE REQUIRED PURPOSE "Votca csg library")
if(DEFINED VOTCA_CSG_VERSION AND NOT VOTCA_CSG_VERSION STREQUAL PROJECT_VERSION)
  message(FATAL_ERROR "Incompatible VOTCA CSG version found (needed ${PROJECT_VERSION}, found ${VOTCA_CSG_VERSION})")
endif()
message(STATUS "Found VOTCA_CSG: ${VOTCA_CSG_DIR}")

if(USE_CUDA)
  find_package(CUDA REQUIRED)
  set_package_properties(CUDA PROPERTIES TYPE REQUIRED PURPOSE "Used for GPU acceleration")
endif(USE_CUDA)

find_package(Kokkos 3.0 CONFIG)
set_package_properties(Kokkos PROPERTIES TYPE OPTIONAL PURPOSE "Used to Ewald solver")

########################################################################
# Checks what linear algebra packages are installed                    #
########################################################################



########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFileCXX)
foreach(HEADER algorithm)
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
add_subdirectory(include/votca/xtp)
add_subdirectory(scripts)
add_subdirectory(src)
add_subdirectory(share)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
endif (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
