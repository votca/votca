cmake_minimum_required(VERSION 3.13)

project(votca LANGUAGES CXX)

set(PROJECT_VERSION "2022-dev")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)

include(GNUInstallDirs)
include(FeatureSummary)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

########################################################################
# User input options                                                   #
########################################################################
option(MODULE_BUILD "Build VOTCA modules one-by-one" OFF)
if(NOT MODULE_BUILD AND NOT EXISTS "${PROJECT_SOURCE_DIR}/tools/CMakeLists.txt")
  message(FATAL_ERROR "No CMake files found in tools/ ! Did you forgot to clone with '--recursive'? Run 'git submodule update --init' to fix that.")
endif()

option(BUILD_XTP "Build xtp" OFF)
add_feature_info(BUILD_XTP BUILD_XTP "Build xtp module")
set(VOTCA_SKIP_FIND_PACKAGE VOTCA_TOOLS VOTCA_CSG VOTCA_XTP)

option(ENABLE_TESTING "Build and copy testing stuff" OFF)
add_feature_info(ENABLE_TESTING ENABLE_TESTING "Enable unit tests")
if(ENABLE_TESTING OR ENABLE_REGRESSION_TESTING)
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
  include(CheckCXXCompilerFlag)
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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcheck")
  endif()
endif()

option(INJECT_MARCH_NATIVE "Inject -march=native if compiler supports it" OFF)
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
if(MODULE_BUILD)
  message(STATUS "Doing Module build")
  include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/BuildModules.cmake)
  return()
endif()

option(BUILD_OWN_GROMACS "Build our own gromacs" OFF)
add_feature_info( BUILD_OWN_GROMACS BUILD_OWN_GROMACS "Build an internal version of gromacs")
if(BUILD_OWN_GROMACS)
  include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/BuildGromacs.cmake)
  list(APPEND VOTCA_SKIP_FIND_PACKAGE GROMACS)
endif()

option(BUILD_OWN_LIBINT "Build our own lib" OFF)
add_feature_info(BUILD_OWN_LIBINT BUILD_OWN_LIBINT "Build an internal version of libint")
if(BUILD_OWN_LIBINT)
  include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/BuildLibint.cmake)
endif()

########################################################################
#Find external packages
########################################################################

#neat trick to not call find_package(VOTCA_*)
macro(find_package)
  if(NOT "${ARGV0}" IN_LIST VOTCA_SKIP_FIND_PACKAGE)
    _find_package(${ARGV})
  endif()
endmacro()

if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
  set_package_properties(Git PROPERTIES TYPE OPTIONAL PURPOSE "Generated version for development version")
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)

######################################
# Include the following subdirectory # 
######################################
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
