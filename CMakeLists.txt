cmake_minimum_required(VERSION 3.1)

project(votca)

set(PROJECT_VERSION "1.5_rc1")

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

option(ENABLE_TESTING "Build and copy testing stuff" OFF)
if(ENABLE_TESTING)
  enable_testing()
endif()

option(BUILD_SHARED_LIBS "Build shared libs" ON)

option(ENABLE_COVERAGE_BUILD "Do a coverage build" OFF)
if(ENABLE_COVERAGE_BUILD)
    message(STATUS "Enabling coverage build")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage -O0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -O0")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()

option(WITH_SQLITE3 "Use SQLite3, at this point only needed for charge transport parts!" ON)
if(NOT WITH_SQLITE3 AND BUILD_XTP )
  message(FATAL_ERROR "xtp needs tools to be build with sqlite3")
endif()

option(ENABLE_RPATH_INJECT "Inject link and install libdir into executables" OFF)
if(ENABLE_RPATH_INJECT)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
  set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_FULL_LIBDIR})
endif(ENABLE_RPATH_INJECT)

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
endif()

########################################################################
#Find external packages
########################################################################
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)

######################################
# Include the following subdirectory # 
######################################
add_subdirectory(tools)
set(ENABLED_VOTCA_PACKAGES tools)
set(VOTCA_TOOLS_LIBRARY votca_tools CACHE STRING "tools library name")
set(VOTCA_TOOLS_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/tools/include;${CMAKE_CURRENT_BINARY_DIR}/tools/include"
  CACHE PATH "tools include path")
set(VOTCA_COMPARE ${CMAKE_BINARY_DIR}/tools/scripts/votca_compare)
add_subdirectory(csg)
list(APPEND ENABLED_VOTCA_PACKAGES csg)
set(VOTCA_SHARE ${CMAKE_SOURCE_DIR}/csg/share)
set(CSG_INVERSE ${CMAKE_BINARY_DIR}/csg/scripts/csg_inverse)
set(CSG_CALL ${CMAKE_BINARY_DIR}/csg/scripts/csg_call)
set(VOTCA_PATH "${CMAKE_BINARY_DIR}/csg/scripts:${CMAKE_BINARY_DIR}/csg/src/tools")
if(BUILD_OWN_GROMACS)
  add_dependencies(votca_csg Gromacs_build)
  set(VOTCA_PATH "${VOTCA_PATH}:${GROMACS_PATH}")
endif()
set(VOTCA_CSG_DEFAULTS ${CMAKE_BINARY_DIR}/csg/share/xml/csg_defaults.xml)
add_subdirectory(csg-tutorials)
list(APPEND ENABLED_VOTCA_PACKAGES csg-tutorials)
if(BUILD_CSG_MANUAL)
  add_subdirectory(csg-manual)
endif()
set(VOTCA_CSG_LIBRARY votca_csg CACHE STRING "csg library name")
set(VOTCA_CSG_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/csg/include
	CACHE PATH "csg include path")
if(BUILD_CSGAPPS)
  add_subdirectory(csgapps)
  list(APPEND ENABLED_VOTCA_PACKAGES csgapps)
endif()
if(BUILD_XTP)
  add_subdirectory(xtp)
  list(APPEND ENABLED_VOTCA_PACKAGES xtp)
endif()

find_package(Doxygen)
if (DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
  add_subdirectory(share/doc)
endif()
add_subdirectory(share/format)

include(FeatureSummary)
feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
