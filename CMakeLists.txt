cmake_minimum_required(VERSION 2.6)

project(votca)

set(PROJECT_VERSION "1.3-dev")

########################################################################
# User input options                                                   #
########################################################################
option(BUILD_CSGAPPS "Build csgapps" OFF)

######################################
# Include the following subdirectory # 
######################################
add_subdirectory(tools)
set(VOTCA_TOOLS_LIBRARY votca_tools)
set(VOTCA_TOOLS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tools/include)
add_subdirectory(csg)
set(VOTCA_CSG_LIBRARY votca_csg)
set(VOTCA_CSG_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/csg/include)
if(BUILD_CSGAPPS)
  add_subdirectory(csgapps)
endif()
