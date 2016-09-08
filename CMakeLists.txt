cmake_minimum_required(VERSION 2.8.3)

project(votca-tools)

set(PROJECT_VERSION "1.4-dev")
string(REGEX REPLACE "^[1-9]+\\.([1-9]+).*$" "\\1" SOVERSION "${PROJECT_VERSION}")
if (NOT ${SOVERSION} MATCHES "[1-9]+")
  message(FATAL_ERROR "Could not determind SOVERSION from ${PROJECT_VERSION}")
endif (NOT ${SOVERSION} MATCHES "[1-9]+")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
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

#this has to be the first include directive
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

########################################################################
#Find external packages
########################################################################
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)

find_package(Threads REQUIRED)
set(THREAD_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

find_package(Boost 1.39.0 REQUIRED COMPONENTS program_options filesystem system )
include_directories(${Boost_INCLUDE_DIRS})
set (BOOST_CFLAGS_PKG "-I${Boost_INCLUDE_DIRS}")
set(BOOST_LIBS_PKG "-L${Boost_LIBRARY_DIRS}")
foreach(_blib ${Boost_LIBRARIES})
  string(REGEX REPLACE ".*/lib([^/]*)\\.[^.]*$" "-l\\1" _blib ${_blib})
  set(BOOST_LIBS_PKG "${BOOST_LIBS_PKG} ${_blib}")
endforeach(_blib)

option(BUILD_MANPAGES "Build manpages (might lead to problem on system without rpath" ON)
#define this target here, so that individual man pages can append to it.
add_custom_target(manpages ALL)

########################################################################
# Checks what linear algebra packages are installed                    #
########################################################################
find_package(GSL)
find_package(MKL)

if (MKL_FOUND AND CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  message("-- MKL found, overloading uBLAS prod()")
  set(MKL TRUE) #used in votca_tools_config.h
elseif(GSL_FOUND)
  message("-- GSL found, overloading uBLAS prod()")
  include_directories(${GSL_INCLUDE_DIRS})
  set(GSL_PKG "gsl")
  set(GSL TRUE) #used in votca_xtp_config.h
endif()

if(MKL_FOUND)
  set(MKL_URL "http://software.intel.com/sites/default/files/managed/76/8c/intel-mkl-and-boost-example.zip")
  set(MKL_MD5 B957EE63A2167AD21B650B47726CD55C )
  set(MKL_FILE "intel-mkl-and-boost-example.zip")
  
  include(ExternalProject)
  ExternalProject_add(intel
          PREFIX intel
          URL "${MKL_URL}" URL_MD5 "${MKL_MD5}"
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          INSTALL_COMMAND ""
          )
  set(MKL_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/intel/src/intel/MKL_and_boost_example/headers")
  include_directories(${MKL_INCLUDE_DIR})
  file(GLOB_RECURSE MKL_HEADERS ${MKL_INCLUDE_DIR}/*.hpp )
  install(FILES ${MKL_HEADERS} DESTINATION include/votca/tools)
endif(MKL_FOUND)

option(WITH_FFTW "Use FFTW3, disabling leads to reduced functionality!" ON)
if (WITH_FFTW)
  find_package(FFTW3)
  if(NOT FFTW3_FOUND)
    message(FATAL_ERROR "FFTW3 not found,  make sure you have also installed the fftw3 and it's dev package (it can be disable with -DWITH_FFTW=OFF)")
  endif(NOT FFTW3_FOUND)
  include_directories(${FFTW3_INCLUDE_DIRS})
  set(FFTW3_PKG "fftw3")
else(WITH_FFTW)
  #used in votca_config.h
  set(NOFFTW TRUE)
  set(FFTW3_PKG)
endif(WITH_FFTW)

#configure_file(config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

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
add_subdirectory(src)
add_subdirectory(include/votca/tools)
add_subdirectory(scripts)
add_subdirectory(share/man)
add_subdirectory(share/doc)

configure_file(${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
