cmake_minimum_required(VERSION 2.6)

project(votca-csg)

set(PROJECT_VERSION "1.2.0-dev")
set(PROJECT_CONTACT "bugs@votca.org")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

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

########################################################################
#Find external packages
########################################################################
find_package(Mercurial)
if (MERCURIAL_FOUND)
  MERCURIAL_HG_INFO(${PROJECT_SOURCE_DIR} CSG)
  MESSAGE("Current revision is ${CSG_HG_ID}")
  set (CSG_HG_ID "hgid: ${CSG_HG_ID}")
else(MERCURIAL_FOUND)
  set (CSG_HG_ID)
endif(MERCURIAL_FOUND)

find_package(TXT2TAGS)
find_package(VOTCA_TOOLS REQUIRED)
include_directories(${VOTCA_TOOLS_INCLUDE_DIRS})
find_package(UnixCommands)

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
foreach(HEADER algorithm fstream iomanip iostream list map numeric sstream stdexcept string vector)
  check_include_file_cxx(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)
include_directories(${CMAKE_SOURCE_DIR}/include)

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
add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(scripts)
add_subdirectory(share)

if (TXT2TAGS_FOUND)
  add_custom_command(OUTPUT CHANGELOG
    COMMAND ${TXT2TAGS_EXECUTABLE} -q -t txt -i ${CMAKE_CURRENT_SOURCE_DIR}/share/CHANGELOG.t2t -o CHANGELOG
    COMMENT "Building CHANGELOG using txt2tags")
else (TXT2TAGS_FOUND)
  add_custom_command(OUTPUT CHANGELOG
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/share/CHANGELOG.t2t CHANGELOG 
    COMMENT "Building CHANGELOG using copy")
endif (TXT2TAGS_FOUND)
add_custom_target(changelog ALL DEPENDS CHANGELOG)

configure_file(${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
