cmake_minimum_required(VERSION 3.10)

project(votca)

set(PROJECT_VERSION "1.6_rc1")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)

include(GNUInstallDirs)
########################################################################
# User input options                                                   #
########################################################################
option(BUILD_CSGAPPS "Build csgapps" OFF)
option(BUILD_XTP "Build xtp" OFF)
option(BUILD_CSG_MANUAL "Build csg-manual" OFF)
set(VOTCA_SKIP_FIND_PACKAGE VOTCA_TOOLS VOTCA_CSG)

option(ENABLE_TESTING "Build and copy testing stuff" OFF)
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
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-conversion")
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcheck")
  endif()
  option(ENABLE_WERROR "Inject -Werror" OFF)
  if(ENABLE_WERROR)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  endif()
endif()

########################################################################
#Find external packages
########################################################################
option(MODULE_BUILD "Build VOTCA modules one-by-one" OFF)
if(MODULE_BUILD)
  message(STATUS "Doing Module build")
  include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/BuildModules.cmake)
  return()
endif()

option(BUILD_OWN_GROMACS "Build our own gromacs" OFF)
if(BUILD_OWN_GROMACS)
  include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/BuildGromacs.cmake)
  list(APPEND VOTCA_SKIP_FIND_PACKAGE GROMACS)
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
if(BUILD_CSG_MANUAL)
  add_subdirectory(csg-manual)
endif()
if(BUILD_CSGAPPS)
  add_subdirectory(csgapps)
  list(APPEND ENABLED_VOTCA_PACKAGES csgapps)
endif()
if(BUILD_XTP)
  add_subdirectory(xtp)
  list(APPEND ENABLED_VOTCA_PACKAGES xtp)
  add_subdirectory(xtp-tutorials)
  list(APPEND ENABLED_VOTCA_PACKAGES xtp-tutorials)
endif()

find_package(Doxygen)
if (DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
  add_subdirectory(share/doc)
endif()
add_subdirectory(share/format)

include(FeatureSummary)
feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
