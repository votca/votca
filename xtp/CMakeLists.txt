cmake_minimum_required(VERSION 3.13)
project(votca-xtp LANGUAGES CXX)

set(PROJECT_VERSION "2023-dev")

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(CheckCXXCompilerFlag)

########################################################################
# User input options                                                   #
########################################################################

find_package(Libxc QUIET CONFIG)
set_package_properties(Libxc PROPERTIES TYPE REQUIRED PURPOSE "Provides exchange-correclation functionals")
if(NOT Libxc_FOUND)
  find_package(Libxc MODULE REQUIRED)
endif()
if(NOT TARGET Libint2::int2)
  find_package(libint2 QUIET CONFIG)
  set_package_properties(libint2 PROPERTIES TYPE REQUIRED PURPOSE "Calculates Gaussian integrals")
  if(NOT libint2_FOUND)
    find_package(libint2 2.6 MODULE REQUIRED)
  endif()
endif()
find_package(ecpint REQUIRED)
set_package_properties(ecpint PROPERTIES TYPE REQUIRED PURPOSE "Calculates Gaussian integrals over pseudo potentials")
message(STATUS "Found libecpint: ${ecpint_DIR}")

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
find_package(Python 3 COMPONENTS Interpreter)
set_package_properties(Python PROPERTIES TYPE OPTIONAL PURPOSE "Used to build documentation files")
option(BUILD_MANPAGES "Build manpages" ${Python_FOUND})
add_feature_info(BUILD_MANPAGES BUILD_MANPAGES "Build manpages (disable for cross-compile)")
#define this target here, so that individual man pages can append to it.
if(NOT TARGET manpages)
  add_custom_target(manpages ALL)
endif()

if(USE_CUDA)
  find_package(CUDA REQUIRED)
  set_package_properties(CUDA PROPERTIES TYPE REQUIRED PURPOSE "Used for GPU acceleration")
endif(USE_CUDA)

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
if(VOTCA_SPHINX_DIR)
  add_custom_target(doc-xtp)
  file(MAKE_DIRECTORY ${VOTCA_SPHINX_DIR}/xtp)
endif()
add_subdirectory(include/votca/xtp)
add_subdirectory(scripts)
add_subdirectory(src)
add_subdirectory(share)
