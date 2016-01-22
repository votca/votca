cmake_minimum_required(VERSION 2.6)
project(votca-xtp)

set(PROJECT_VERSION "1.3-dev")
set(PROJECT_CONTACT "bugs@votca.org")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)


enable_language(CXX)

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(CheckCXXCompilerFlag)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


########################################################################
# User input options                                                   #
########################################################################


find_package(LIBXC)
if (LIBXC_FOUND)
	include_directories(${LIBXC_INCLUDE_DIRS})
	set (LIBXC TRUE )
else()
	set(LIBXC_LIBRARIES)
endif()


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
#option(WITH_DIPRO "Build mergeorbitals and DIPRO, which use fortran" ON)

#for votca_xtp_config.h
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
########################################################################
#Find external packages
########################################################################
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.hg)
  find_package(Mercurial)
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.hg)

if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)

find_package(Boost 1.48.0 REQUIRED COMPONENTS program_options serialization filesystem timer)
include_directories(${Boost_INCLUDE_DIRS})
set (BOOST_LIBRARIES ${Boost_THREAD_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_TIMER_LIBRARY})
set (BOOST_CFLAGS_PKG "-I${Boost_INCLUDE_DIRS}")
string(REGEX REPLACE "^(.*)/lib([^/]*)\\.[^.]*$" "-L\\1 -l\\2" BOOST_LIBS_PKG ${BOOST_LIBRARIES})

option(BUILD_MANPAGES "Build manpages (might lead to problem on system without rpath" OFF)
#define this target here, so that individual man pages can append to it.
if(NOT TARGET manpages)
  add_custom_target(manpages ALL)
endif()




#message(STATUS "BOOST_LIBRARIES='${BOOST_LIBRARIES}'")

find_package(VOTCA_TOOLS REQUIRED)
include_directories(${VOTCA_TOOLS_INCLUDE_DIRS})
find_package(VOTCA_CSG REQUIRED)
include_directories(${VOTCA_CSG_INCLUDE_DIRS})

########################################################################
# Checks what linear algebra packages are installed                    #
########################################################################
option(WITH_LINALG "Checking which BLAS libraries are installed" ON)
if (WITH_LINALG)

  find_package(Eigen)


  if (EIGEN_FOUND) 
    message("-- EIGEN found (not used at the moment)")
    set(EIGEN TRUE) #used in votca_xtp_config.h
    include_directories(${EIGEN_INCLUDE_DIRS})
  else()
    message("-- EIGEN not found")
    set(EIGEN FALSE) #used in votca_xtp_config.h
  endif()

  set(MKL FALSE)
  set(GSL FALSE)

  #find_package(GSL)
  find_package(MKL)


  if (MKL_FOUND AND CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message("-- MKL found, overloading uBLAS prod()")
    set(MKL TRUE) #used in votca_xtp_config.h
  elseif(GSL_FOUND)
    message("-- GSL found, overloading uBLAS prod()")
    include_directories(${GSL_INCLUDE_DIRS})
    set(GSL_PKG "gsl")
    set(GSL TRUE) #used in votca_xtp_config.h
  endif()

endif(WITH_LINALG)

#configure_file(config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

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
add_subdirectory(include/votca/xtp)
add_subdirectory(scripts)
add_subdirectory(src/libxtp)
add_subdirectory(src/tools)
#add_subdirectory(src/dipro)
#add_subdirectory(share/scripts/dipro) 
add_subdirectory(share)  
