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

########################################################################
#Find external packages
########################################################################
find_package(Mercurial)
find_package(Threads REQUIRED)
set(THREAD_LIBRARY ${CMAKE_THREAD_LIBS_INIT})

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFiles)
check_include_files("sys/time.h;sys/resource.h;unistd.h" CMAKE_HEADERS)

include_directories(${CMAKE_SOURCE_DIR}/include)

######################################
# Include the following subdirectory # 
######################################
add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(scripts)

