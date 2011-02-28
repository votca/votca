cmake_minimum_required(VERSION 2.6)

project(votca-tools)

set(PROJECT_VERSION "1.2.0-dev")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)


enable_language(CXX)


######################################################################
# cpack seetings
#####################################################################
#SET(CPACK_PACKAGE_VERSION_MAJOR "1")
#SET(CPACK_PACKAGE_VERSION_MINOR "2")
#SET(CPACK_PACKAGE_VERSION_PATCH "0")
#set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
#set(CPACK_GENERATOR "")
#set(CPACK_SOURCE_GENERATOR "TGZ")
#set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
#set(CPACK_PACKAGE_VENDOR "espresso-pp.de")
#set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PROJECT_NAME}-${PROJECT_VERSION}")
#set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ESPResSo++ is a highly versatile software package for the scientific simulation ")
#set(CPACK_SOURCE_IGNORE_FILES ".*~;/\\\\.hg/;\\\\.hgignore;/CMakeFiles/;CMakeCache\\\\.txt$;cmake_install\\\\.cmake$;CPack.*\\\\.cmake$;Makefile$;ESPRC$")
#if (EXTERNAL_BOOST)
#  set(CPACK_SOURCE_IGNORE_FILES "${CPACK_SOURCE_IGNORE_FILES};contrib/boost/")
#endif (EXTERNAL_BOOST)
#include(CPack)


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
find_package(TXT2TAGS)

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFile)
foreach(HEADER assert.h math.h stdio.h stdlib.h string.h)
  check_include_file(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

include(CheckIncludeFileCXX)
foreach(HEADER cmath fstream functional iostream limits list map ostream sstream stack stdexcept string vector)
  check_include_file_cxx(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)
include_directories(${CMAKE_SOURCE_DIR}/include)

include(CheckLibraryExists)
foreach(FUNC sqrt)
  check_library_exists(m ${FUNC} "" FOUND_${FUNC})
  if(NOT FOUND_${FUNC})
    message(FATAL_ERROR "Could not find needed math function - ${FUNC}")
  endif(NOT FOUND_${FUNC})
endforeach(FUNC)
set(MATH_LIBRARIES m)

######################################
# Include the following subdirectory # 
######################################
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
