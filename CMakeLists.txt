cmake_minimum_required(VERSION 2.6)

project(votca)

#cmake policy CMP0048
#set(PROJECT_VERSION "1.5-dev")

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

######################################
# Include the following subdirectory # 
######################################
add_subdirectory(tools)
set(VOTCA_TOOLS_LIBRARY votca_tools CACHE STRING "tools library name")
set(VOTCA_TOOLS_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/tools/include;${CMAKE_CURRENT_BINARY_DIR}/tools/include"
  CACHE PATH "tools include path")
add_subdirectory(csg)
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
endif()
if(BUILD_CTP)
  add_subdirectory(kmc)
  add_subdirectory(moo)
  set(VOTCA_MOO_LIBRARY votca_moo CACHE STRING "moo library name")
  set(VOTCA_MOO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/moo/include
	  CACHE PATH "moo include path")
  add_subdirectory(ctp)
  set(VOTCA_CTP_LIBRARY votca_ctp CACHE STRING "ctp library name")
  set(VOTCA_CTP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ctp/include
	  CACHE PATH "ctp include path")
endif()
if(BUILD_XTP)
  add_subdirectory(xtp)
endif()

include(FeatureSummary)
feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
