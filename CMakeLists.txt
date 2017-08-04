cmake_minimum_required(VERSION 3.1)

project(votca)

#cmake policy CMP0048
#set(PROJECT_VERSION "1.5-dev")
set(VOTCA_VERSION "1.5-dev")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)

########################################################################
# User input options                                                   #
########################################################################
option(BUILD_CSGAPPS "Build csgapps" OFF)
option(BUILD_CTP "Build ctp" OFF)
option(BUILD_XTP "Build xtp" OFF)
if(BUILD_XTP AND NOT BUILD_CTP)
  message(FATAL_ERROR "xtp needs ctp to be build")
endif(BUILD_XTP AND NOT BUILD_CTP)
option(BUILD_CSG_MANUAL "Build csg-manual" OFF)

option(ENABLE_TESTING "Build and copy testing stuff" OFF)
if(ENABLE_TESTING)
  enable_testing()
endif()

option(ENABLE_COVERAGE_BUILD "Do a coverage build" OFF)
if(ENABLE_COVERAGE_BUILD)
    message(STATUS "Enabling coverage build")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage -O0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -O0")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()

######################################
# Include the following subdirectory # 
######################################
add_subdirectory(tools)
set(DOXYGEN_PACKAGES tools)
set(VOTCA_TOOLS_LIBRARY votca_tools CACHE STRING "tools library name")
set(VOTCA_TOOLS_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/tools/include;${CMAKE_CURRENT_BINARY_DIR}/tools/include"
  CACHE PATH "tools include path")
add_subdirectory(csg)
list(APPEND DOXYGEN_PACKAGES csg)
set(VOTCA_SHARE ${CMAKE_SOURCE_DIR}/csg/share)
set(VOTCA_CSG_INVERSE ${CMAKE_BINARY_DIR}/csg/scripts/csg_inverse)
set(CSG_CALL ${CMAKE_BINARY_DIR}/csg/scripts/csg_call)
set(CSG_INVERSE ${VOTCA_CSG_INVERSE})
set(VOTCA_PATH "${CMAKE_BINARY_DIR}/csg/scripts:${CMAKE_BINARY_DIR}/csg/src/tools")
add_subdirectory(csg-tutorials)
if(BUILD_CSG_MANUAL)
  add_subdirectory(csg-manual)
endif()
set(VOTCA_CSG_LIBRARY votca_csg CACHE STRING "csg library name")
set(VOTCA_CSG_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/csg/include
	CACHE PATH "csg include path")
if(BUILD_CSGAPPS)
  add_subdirectory(csgapps)
  list(APPEND DOXYGEN_PACKAGES csgapps)
endif()
if(BUILD_CTP)
  add_subdirectory(kmc)
  list(APPEND DOXYGEN_PACKAGES kmc)
  add_subdirectory(moo)
  list(APPEND DOXYGEN_PACKAGES moo)
  set(VOTCA_MOO_LIBRARY votca_moo CACHE STRING "moo library name")
  set(VOTCA_MOO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/moo/include
	  CACHE PATH "moo include path")
  add_subdirectory(ctp)
  list(APPEND DOXYGEN_PACKAGES ctp)
  set(VOTCA_CTP_LIBRARY votca_ctp CACHE STRING "ctp library name")
  set(VOTCA_CTP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ctp/include
	  CACHE PATH "ctp include path")
endif()
if(BUILD_XTP)
  add_subdirectory(xtp)
  list(APPEND DOXYGEN_PACKAGES ctp)
endif()

find_package(Doxygen)
if (DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
  add_subdirectory(share/doc)
endif()

include(FeatureSummary)
feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
