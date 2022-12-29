cmake_minimum_required(VERSION 3.13)

project(votca LANGUAGES CXX)

set(PROJECT_VERSION "2023-dev")
string(REGEX REPLACE "[-.].*$" "" SOVERSION "${PROJECT_VERSION}")
if (NOT ${SOVERSION} MATCHES "^[0-9]+$")
  message(FATAL_ERROR "Could not determine SOVERSION (${SOVERSION}) from ${PROJECT_VERSION}")
endif (NOT ${SOVERSION} MATCHES "^[0-9]+$")
set(PROJECT_CONTACT "https://github.com/votca/votca/issues")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)

include(GNUInstallDirs)
include(FeatureSummary)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

########################################################################
# User input options                                                   #
########################################################################
option(BUILD_XTP "Build xtp" OFF)
add_feature_info(BUILD_XTP BUILD_XTP "Build xtp module")

option(ENABLE_TESTING "Build and copy testing stuff" OFF)
add_feature_info(ENABLE_TESTING ENABLE_TESTING "Enable unit tests")
if(ENABLE_TESTING OR ENABLE_REGRESSION_TESTING)
  find_package(Boost 1.71.0 REQUIRED COMPONENTS unit_test_framework)
  enable_testing()
endif()

option(BUILD_SHARED_LIBS "Build shared libs" ON)

option(ENABLE_COVERAGE_BUILD "Do a coverage build" OFF)
if(ENABLE_COVERAGE_BUILD)
    message(STATUS "Enabling coverage build")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()

option(ENABLE_RPATH_INJECT "Inject link and install libdir into executables" OFF)
add_feature_info(ENABLE_RPATH_INJECT ENABLE_RPATH_INJECT "Inject rpath into executables")
if(ENABLE_RPATH_INJECT)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
  set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_FULL_LIBDIR})
endif(ENABLE_RPATH_INJECT)

option(ENABLE_WARNING_FLAGS "Inject more warning flags" ON)
if(ENABLE_WARNING_FLAGS)
  foreach(_FLAG -Wall -Wextra -Wpedantic -Wshadow -Wconversion) 
    check_cxx_compiler_flag("${_FLAG}" COMPILER_SUPPORTS${_FLAG})
    if(COMPILER_SUPPORTS${_FLAG})
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_FLAG}")
    endif()
  endforeach()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-conversion")
  endif()
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel") #only icpc, not icpx
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcheck -Wno-sign-conversion")
  endif()
  # allow adding extra warning flags at the very end
  if(VOTCA_EXTRA_WARNING_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VOTCA_EXTRA_WARNING_FLAGS}")
  endif()
endif()

option(INJECT_MARCH_NATIVE "Inject -march=native if compiler supports it" ON)
add_feature_info(INJECT_MARCH_NATIVE INJECT_MARCH_NATIVE "Add -march=native to CXX flags")
if(INJECT_MARCH_NATIVE)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_march_native)
  if(COMPILER_SUPPORTS_march_native)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif()
endif()

option(ENABLE_WERROR "Inject -Werror" OFF)
if(ENABLE_WERROR)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()

if(CMAKE_CXX_COMPILER MATCHES "dpcpp")
  message(WARNING "Using the intel dpcpp compiler needs a cutting-edge version of eigen hence dpcpp is currently not offically supported, use icpx instead") 
endif()

set(VOTCA_SPHINX_DIR "${PROJECT_BINARY_DIR}/sphinx" CACHE PATH "Path for Sphinx sources")
file(MAKE_DIRECTORY ${VOTCA_SPHINX_DIR})
set(VOTCA_SPHINX_OUTPUT_DIR "${PROJECT_BINARY_DIR}/sphinx.html" CACHE PATH "Path for Sphinx output")
set(VOTCA_XML_PARSER "${PROJECT_SOURCE_DIR}//share/doc/extract_xml_metadata.py")

########################################################################
#Find external packages
########################################################################
find_package(Boost 1.71.0 REQUIRED
  COMPONENTS program_options filesystem system regex timer)
set_package_properties(Boost PROPERTIES TYPE REQUIRED PURPOSE "Extended C++ libraries")

find_package(Eigen3 3.3.0 NO_MODULE REQUIRED)
set_package_properties(Eigen3 PROPERTIES TYPE REQUIRED PURPOSE "C++ vector data structures")
message(STATUS "Found Eigen3: ${Eigen3_DIR}")

find_package(Threads REQUIRED)
set_package_properties(Threads PROPERTIES TYPE REQUIRED PURPOSE "Used for thread parallelization")

find_package(GROMACS 2019.6 CONFIG NAMES gromacs gromacs_d)
set_package_properties(GROMACS PROPERTIES TYPE RECOMMENDED PURPOSE "Used to read/write gromacs data files")
if(GROMACS_FOUND AND DEFINED GROMACS_VERSION AND GROMACS_VERSION VERSION_GREATER_EQUAL "2020")
  message(WARNING "Gromacs-2020 and above have no support for tabulated interactions, that are needed for coarse-graining (see and comment on https://gitlab.com/gromacs/gromacs/-/issues/1347)")
endif()
# prior gmx-2021, their exported target had some bugs, that we fix below
# - missing GMX_DOUBLE define - only gets exported as GROMACS_DEFINITIONS
# - missing namespace
if(GROMACS_FOUND AND DEFINED GROMACS_VERSION AND GROMACS_VERSION VERSION_LESS "2021")
  set_target_properties(libgromacs PROPERTIES INTERFACE_COMPILE_OPTIONS ${GROMACS_DEFINITIONS})
  add_library(Gromacs::libgromacs ALIAS libgromacs)
endif()

# iie methods needs python 3.5+
find_package(Python 3.5 REQUIRED COMPONENTS Interpreter)
set_package_properties(Python PROPERTIES TYPE OPTIONAL PURPOSE "Used in csg scripts and docu")
option(BUILD_MANPAGES "Build manpages (might lead to problem on system without rpath" ${Python_FOUND})
add_feature_info(BUILD_MANPAGES BUILD_MANPAGES "Build manpages (disable for cross-compile)")

set(MATH_LIBRARIES "m" CACHE STRING "math library")
mark_as_advanced( MATH_LIBRARIES )

if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
  set_package_properties(Git PROPERTIES TYPE OPTIONAL PURPOSE "Generated version for development version")
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)

option(BUILD_OWN_GROMACS "Build our own gromacs" OFF)
add_feature_info( BUILD_OWN_GROMACS BUILD_OWN_GROMACS "Build an internal version of gromacs")
if(BUILD_OWN_GROMACS)
  include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/BuildGromacs.cmake)
  #neat trick to not call find_package(GROMACS)
  macro(find_package)
    if(NOT "${ARGV0}" STREQUAL GROMACS)
      _find_package(${ARGV})
    endif()
  endmacro()
endif()

option(BUILD_OWN_LIBINT "Build our own lib" OFF)
add_feature_info(BUILD_OWN_LIBINT BUILD_OWN_LIBINT "Build an internal version of libint")
if(BUILD_OWN_LIBINT)
  include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/BuildLibint.cmake)
endif()
######################################
# Include the following subdirectory # 
######################################
configure_file(${CMAKE_SOURCE_DIR}/CMakeModules/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

if(BUILD_MANPAGES)
  #define this target here, so that individual man pages can append to it.
  add_custom_target(manpages ALL)
endif()

add_subdirectory(tools)
set(ENABLED_VOTCA_PACKAGES tools)
add_subdirectory(csg)
list(APPEND ENABLED_VOTCA_PACKAGES csg)
add_subdirectory(csg-tutorials)
list(APPEND ENABLED_VOTCA_PACKAGES csg-tutorials)
if(BUILD_XTP)
  add_subdirectory(xtp)
  list(APPEND ENABLED_VOTCA_PACKAGES xtp)
  add_subdirectory(xtp-tutorials)
  list(APPEND ENABLED_VOTCA_PACKAGES xtp-tutorials)
endif()

add_subdirectory(share/doc)
add_subdirectory(share/format)

feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
