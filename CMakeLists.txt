cmake_minimum_required(VERSION 3.10)

project(votca-csg)

set(PROJECT_VERSION "1.6-dev")
set(PROJECT_CONTACT "bugs@votca.org")
string(REGEX REPLACE "^[1-9]+\\.([1-9]+).*$" "\\1" SOVERSION "${PROJECT_VERSION}")
if (NOT ${SOVERSION} MATCHES "[1-9]+")
  message(FATAL_ERROR "Could not determind SOVERSION from ${PROJECT_VERSION}")
endif (NOT ${SOVERSION} MATCHES "[1-9]+")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_definitions(-DDEBUG)
endif(CMAKE_BUILD_TYPE STREQUAL Debug)

enable_language(CXX)

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

option(ENABLE_TESTING "Build and enable testing stuff" OFF)
if(ENABLE_TESTING)
  enable_testing()
endif(ENABLE_TESTING)

########################################################################
#Find external packages
########################################################################
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)

find_package(Threads REQUIRED)
option(BUILD_MANPAGES "Build manpages" ON)
if (BUILD_MANPAGES)
  find_package(TXT2TAGS)
  find_package(UnixCommands)
else (BUILD_MANPAGES)
  #we use TXT2TAGS_FOUND conditionals in other CMakeLists.txt
  #so that TXT2TAGS is never required
  set(TXT2TAGS_FOUND)
endif (BUILD_MANPAGES)

find_package(Boost 1.53.0 REQUIRED COMPONENTS program_options filesystem system )
find_package(Eigen3 3.3.0 NO_MODULE REQUIRED)
message(STATUS "Found Eigen3: ${Eigen3_DIR}")

find_package(VOTCA_TOOLS NO_MODULE REQUIRED)
message(STATUS "Found VOTCA_TOOLS: ${VOTCA_TOOLS_DIR}")
if(DEFINED VOTCA_TOOLS_VERSION AND NOT VOTCA_TOOLS_VERSION STREQUAL PROJECT_VERSION)
  message(FATAL_ERROR "Incompatible VOTCA Tools version found (needed ${PROJECT_VERSION}, found ${VOTCA_TOOLS_VERSION})")
endif()

option(WITH_GMX "Build gromacs reader/writer, disabling leads to reduced functionality!" ON)
if (WITH_GMX)
  find_package(GROMACS 20160000 REQUIRED)
  if(DEFINED GROMACS_VERSION AND GROMACS_VERSION VERSION_GREATER_EQUAL "20200000")
    message(WARNING "Gromacs-2020 and above have no support for tabulated interactions, that are needed for coarse-graining (see and comment on https://redmine.gromacs.org/issues/1347)")
  endif()
endif(WITH_GMX) 

if(ENABLE_TESTING)
  set(REGRESSIONTEST_TOLERANCE "5e-5" CACHE STRING "Tolerance for the regression tests")
  mark_as_advanced(REGRESSIONTEST_TOLERANCE)
  find_package(UnixCommands)
endif(ENABLE_TESTING)

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFile)
foreach(HEADER assert.h)
  check_include_file(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

include(CheckIncludeFileCXX)
foreach(HEADER algorithm fstream iomanip iostream list map numeric sstream stdexcept string vector cstdlib)
  check_include_file_cxx(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

set(MATH_LIBRARIES "m" CACHE STRING "math library")
mark_as_advanced( MATH_LIBRARIES )
include(CheckLibraryExists)
foreach(FUNC sqrt)
  check_library_exists(${MATH_LIBRARIES} ${FUNC} "" FOUND_${FUNC}_${MATH_LIBRARIES})
  if(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
    message(FATAL_ERROR "Could not find needed math function - ${FUNC}")
  endif(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
endforeach(FUNC)

######################################
# Include the following subdirectory # 
######################################
if(NOT TARGET manpages)
  add_custom_target(manpages ALL)
endif()
add_subdirectory(src)
add_subdirectory(scripts)
add_subdirectory(include/votca/csg)
add_subdirectory(share)

# needs to happen after src and scripts subdir
get_property(VOTCA_BINARIES TARGET votca_csg PROPERTY BINARIES)
list(SORT VOTCA_BINARIES)
configure_file(src/libcsg/VOTCA_CSGConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/src/libcsg/VOTCA_CSGConfig.cmake" @ONLY)
# for >=cmake-3.12 use list(JOIN VOTCA_BINARIES " " VOTCA_BINARIES)
string(REPLACE ";" " " VOTCA_BINARIES "${VOTCA_BINARIES}")
configure_file(scripts/csg-completion.bash.in "${CMAKE_CURRENT_BINARY_DIR}/scripts/csg-completion.bash" @ONLY)

configure_file(${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall-csg COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
if(NOT TARGET uninstall)
  add_custom_target(uninstall)
endif()
add_dependencies(uninstall uninstall-csg)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  include(FeatureSummary)
  feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
endif (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
