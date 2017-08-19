cmake_minimum_required(VERSION 2.6)

project(csgapps)

set(PROJECT_VERSION "1.5-dev")

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

#needed by csg
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
	message(FATAL_ERROR "Could not find a C++-11 compiler")
endif()

########################################################################
# User input options                                                   #
########################################################################
option(ENABLE_TESTING "Build and enable testing stuff" OFF)
if(ENABLE_TESTING)
  enable_testing()
endif(ENABLE_TESTING)

########################################################################
#Find external packages
########################################################################

find_package(Boost 1.39.0 REQUIRED COMPONENTS program_options )
include_directories(${Boost_INCLUDE_DIRS})
set (BOOST_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})

find_package(VOTCA_TOOLS REQUIRED)
include_directories(${VOTCA_TOOLS_INCLUDE_DIRS})
find_package(VOTCA_CSG REQUIRED)
include_directories(${VOTCA_CSG_INCLUDE_DIRS})

file(GLOB_RECURSE FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} */*.cc)
foreach(DIR ${FILES})
  string(REGEX MATCH ^[^/]* DIR ${DIR})
  list(APPEND DIRS ${DIR})
endforeach(DIR)
list(REMOVE_DUPLICATES DIRS)
foreach(DIR ${DIRS})
  file(GLOB SRC ${DIR}/*.cc)
  add_executable(${DIR} ${SRC})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DIR}) #needed for out-of-source build
  target_link_libraries(${DIR} ${VOTCA_CSG_LIBRARIES} ${VOTCA_TOOLS_LIBRARIES} ${BOOST_LIBRARIES})
  install(TARGETS ${DIR} RUNTIME DESTINATION bin)
  set_target_properties(${DIR} PROPERTIES OUTPUT_NAME ${DIR}/csg_${DIR})
  if(ENABLE_TESTING)
    add_test(${DIR}Help ${DIR}/csg_${DIR} --help)
    # run tests for tools and csg as well for coverage
    set_tests_properties(${PROG}Help PROPERTIES LABELS "csgapps;csg;tools;votca")
  endif(ENABLE_TESTING)
endforeach(DIR)

configure_file(${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall-csgapps COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
if(NOT TARGET uninstall)
  add_custom_target(uninstall)
endif()
add_dependencies(uninstall uninstall-csgapps)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  include(FeatureSummary)
  feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
endif (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
