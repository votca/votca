cmake_minimum_required(VERSION 2.6)

project(votca-ctp)

set(PROJECT_VERSION "0.0-dev")
set(PROJECT_CONTACT "bugs@votca.org")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)


enable_language(CXX)

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(CheckCXXCompilerFlag)

########################################################################
# User input options                                                   #
########################################################################
option(BUILD_SHARED_LIBS "Build shared libs" ON)
if (NOT DEFINED LIB)
#  if (CMAKE_SIZEOF_VOID_P MATCHES "8")
#    set( LIB_SUFFIX "64" )
#  endif(CMAKE_SIZEOF_VOID_P MATCHES "8")
#  set(LIB "lib${LIB_SUFFIX}")
  set(LIB "lib")
endif(NOT DEFINED LIB)
if (NOT DEFINED MAN)
  set(MAN "share/man")
endif(NOT DEFINED MAN)
if (NOT DEFINED DATA)
  set(DATA "share/votca")
endif(NOT DEFINED DATA)

include_directories(${CMAKE_SOURCE_DIR}/include)
########################################################################
#Find external packages
########################################################################
option(EXTERNAL_BOOST "Use external boost" ON)
if (EXTERNAL_BOOST)
  find_package(Boost 1.39.0 COMPONENTS program_options )
  if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost not found,  make sure you have also installed boost and it's dev package or use our internal replacement, if you have built tools with -DEXTERNAL_BOOST=OFF and set it here as well")
  endif(NOT Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  set (BOOST_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
  set (BOOST_PKG)
  set (BOOST_CFLAGS_PKG "-I${Boost_INCLUDE_DIRS}")
  string(REGEX REPLACE "^(.*)/lib([^/]*)\\.[^.]*$" "-L\\1 -l\\2" BOOST_LIBS_PKG ${BOOST_LIBRARIES})
else(EXTERNAL_BOOST)
  find_package(VOTCA_BOOST REQUIRED)
  include_directories(${VOTCA_BOOST_INCLUDE_DIRS})
  set (BOOST_LIBRARIES ${VOTCA_BOOST_LIBRARIES})
  set (BOOST_PKG "libvotca_boost")
  set (BOOST_INCLUDE_DIR_PKG)
  set (BOOST_LIB_PKG)
  set (BOOST_LIB_DIR_PKG)
endif(EXTERNAL_BOOST)

find_package(VOTCA_TOOLS REQUIRED)
include_directories(${VOTCA_TOOLS_INCLUDE_DIRS})
find_package(VOTCA_CSG REQUIRED)
include_directories(${VOTCA_CSG_INCLUDE_DIRS})
find_package(VOTCA_MOO REQUIRED)
include_directories(${VOTCA_MOO_INCLUDE_DIRS})


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
foreach(HEADER algorithm)
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
add_subdirectory(include/votca/ctp)
add_subdirectory(scripts)
add_subdirectory(src/libctp)
add_subdirectory(src/tools)
add_subdirectory(src/dipro)
add_subdirectory(share/scripts/dipro) 
add_subdirectory(share)  
