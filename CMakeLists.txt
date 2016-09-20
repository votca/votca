cmake_minimum_required(VERSION 2.6)
project(votca-xtp)

set(PROJECT_VERSION "1.4-dev")
set(PROJECT_CONTACT "bugs@votca.org")
string(REGEX REPLACE "^[1-9]+\\.([1-9]+).*$" "\\1" SOVERSION "${PROJECT_VERSION}")
if (NOT ${SOVERSION} MATCHES "[1-9]+")
  message(FATAL_ERROR "Could not determind SOVERSION from ${PROJECT_VERSION}")
endif (NOT ${SOVERSION} MATCHES "[1-9]+")

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
if (NOT DEFINED XTP_MANUAL_LOCATION)
  set(XTP_MANUAL_LOCATION "share/votca/manual")
endif(NOT DEFINED XTP_MANUAL_LOCATION)
#option(WITH_DIPRO "Build mergeorbitals and DIPRO, which use fortran" ON)

#for votca_config.h
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
########################################################################
#Find external packages
########################################################################
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



set(MKL FALSE)
set(GSL FALSE)

find_package(GSL)
find_package(MKL)


if (MKL_FOUND AND CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
message("-- MKL found, overloading uBLAS prod()")
set(MKL TRUE) #used in votca_config.h
elseif(GSL_FOUND)
message("-- GSL found, overloading uBLAS prod()")
include_directories(${GSL_INCLUDE_DIRS})
set(GSL_PKG "gsl")
set(GSL TRUE) #used in votca_config.h
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
  install(FILES ${MKL_HEADERS} DESTINATION include/votca/xtp)
endif(MKL_FOUND)



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

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/manual/xtp-manual.pdf COMMAND make -C ${CMAKE_SOURCE_DIR}/manual
	XTP_SCRIPTS="${XTP_SCRIPTS}"
	CSG_PROPERTY="${VOTCA_CSG_PROPERTY}"
	VER=${PROJECT_VERSION}
	ROOTDIR="${CMAKE_CURRENT_SOURCE_DIR}/manual"
	XTP_PROGS="${XTP_PROGS}" VERBATIM)
add_custom_target(manual DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/manual/xtp-manual.pdf)

add_custom_target(cleanmanual COMMAND make clean -C ${CMAKE_SOURCE_DIR}/manual
	XTP_SCRIPTS="${XTP_SCRIPTS}"
	CSG_PROPERTY="${VOTCA_CSG_PROPERTY}"
	VER=${PROJECT_VERSION}
	ROOTDIR="${CMAKE_CURRENT_SOURCE_DIR}/manual"
	XTP_PROGS="${XTP_PROGS}" VERBATIM)

add_custom_target(installmanual COMMAND make install -C ${CMAKE_SOURCE_DIR}/manual
	XTP_SCRIPTS="${XTP_SCRIPTS}"
	CSG_PROPERTY="${VOTCA_CSG_PROPERTY}"
	VER=${PROJECT_VERSION}
	ROOTDIR="${CMAKE_CURRENT_SOURCE_DIR}/manual"
	XTP_PROGS="${XTP_PROGS}"
	XTP_MANUAL_LOCATION="${CMAKE_INSTALL_PREFIX}/${XTP_MANUAL_LOCATION}" VERBATIM)
add_dependencies(installmanual manual)
