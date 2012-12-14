cmake_minimum_required(VERSION 2.6)

project(votca-csg)

set(PROJECT_VERSION "1.3-dev")
set(PROJECT_CONTACT "bugs@votca.org")
string(REGEX REPLACE "^[1-9]+\\.([1-9]+).*$" "\\1" SOVERSION "${PROJECT_VERSION}")
if (NOT ${SOVERSION} MATCHES "[1-9]+")
  message(FATAL_ERROR "Could not determind SOVERSION from ${PROJECT_VERSION}")
endif (NOT ${SOVERSION} MATCHES "[1-9]+")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
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

########################################################################
# User input options                                                   #
########################################################################
option(BUILD_SHARED_LIBS "Build shared libs" ON)
if (NOT DEFINED LIB)
  set(LIB "lib")
endif(NOT DEFINED LIB)
if (NOT DEFINED MAN)
  set(MAN "share/man")
endif(NOT DEFINED MAN)
if (NOT DEFINED DATA)
  set(DATA "share/votca")
endif(NOT DEFINED DATA)

# this has to be the include specification!
include_directories(${CMAKE_SOURCE_DIR}/include)

########################################################################
#Find external packages
########################################################################
option(BUILD_MANPAGES "Build manpages" ON)
if (BUILD_MANPAGES)
  find_package(TXT2TAGS)
  find_package(UnixCommands)
else (BUILD_MANPAGES)
  #we use TXT2TAGS_FOUND conditionals in other CMakeLists.txt
  #so that TXT2TAGS is never required
  set(TXT2TAGS_FOUND)
endif (BUILD_MANPAGES)

find_package(Boost 1.39.0 REQUIRED COMPONENTS program_options )
include_directories(${Boost_INCLUDE_DIRS})
set (BOOST_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
set (BOOST_CFLAGS_PKG "-I${Boost_INCLUDE_DIRS}")
string(REGEX REPLACE "^(.*)/lib([^/]*)\\.[^.]*$" "-L\\1 -l\\2" BOOST_LIBS_PKG ${BOOST_LIBRARIES})

find_package(VOTCA_TOOLS REQUIRED)
include_directories(${VOTCA_TOOLS_INCLUDE_DIRS})

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
foreach(HEADER algorithm fstream iomanip iostream list map numeric sstream stdexcept string vector cstdlib)
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
add_custom_target(manpages ALL)
add_subdirectory(scripts)
add_subdirectory(src)
add_subdirectory(include/votca/csg)
add_subdirectory(share)

file(GLOB VOTCA_CMAKE_MODULES CMakeModules/FindVOTCA* CMakeModules/FindSQLITE3.cmake CMakeModules/FindBoost.cmake)
install(FILES ${VOTCA_CMAKE_MODULES} DESTINATION ${DATA}/template/CMakeModules)

configure_file(${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
