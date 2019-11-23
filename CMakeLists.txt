cmake_minimum_required(VERSION 3.10)

project(csgapps)

set(PROJECT_VERSION "1.6-dev")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)

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
include(GNUInstallDirs)

option(ENABLE_TESTING "Build and enable testing stuff" OFF)
if(ENABLE_TESTING)
  enable_testing()
endif(ENABLE_TESTING)

########################################################################
#Find external packages
########################################################################
find_package(VOTCA_CSG REQUIRED NO_MODULE)
if(DEFINED VOTCA_CSG_VERSION AND NOT VOTCA_CSG_VERSION STREQUAL PROJECT_VERSION)
  message(FATAL_ERROR "Incompatible VOTCA CSG version found (needed ${PROJECT_VERSION}, found ${VOTCA_CSG_VERSION})")
endif()
message(STATUS "Found VOTCA_CSG: ${VOTCA_CSG_DIR}")

file(GLOB_RECURSE FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} */*.cc)
foreach(DIR ${FILES})
  string(REGEX MATCH ^[^/]* DIR ${DIR})
  list(APPEND DIRS ${DIR})
endforeach(DIR)
list(REMOVE_DUPLICATES DIRS)
foreach(DIR ${DIRS})
  file(GLOB SRC ${DIR}/*.cc)
  add_executable(${DIR} ${SRC})
  target_link_libraries(${DIR} VOTCA::votca_csg)
  install(TARGETS ${DIR} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  set_target_properties(${DIR} PROPERTIES OUTPUT_NAME csg_${DIR} RUNTIME_OUTPUT_DIRECTORY ${DIR})
  if(ENABLE_TESTING)
    add_test(${DIR}Help ${DIR}/csg_${DIR} --help)
    # run tests for tools and csg as well for coverage
    set_tests_properties(${DIR}Help PROPERTIES LABELS "csgapps;csg;tools;votca")
  endif(ENABLE_TESTING)
endforeach(DIR)

configure_file(${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall-csgapps COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
if(NOT TARGET uninstall)
  add_custom_target(uninstall)
endif()
add_dependencies(uninstall uninstall-csgapps)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  include(FeatureSummary)
  feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
endif (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
