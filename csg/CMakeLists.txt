cmake_minimum_required(VERSION 3.13)

project(votca-csg LANGUAGES CXX)

set(PROJECT_VERSION "2023-dev")
########################################################################
#Find external packages
########################################################################
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
