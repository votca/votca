# - Finds parts of gromacs
# Find the native gromacs compents headers and libraries.
#
#  GROMACS_INCLUDE_DIRS - where to find gromacs headers.
#  GROMACS_LIBRARIES    - List of libraries when used by gromacs.
#  GROMACS_FOUND        - True if all gromacs componets were found.
#  GROMACS_DEFINITIONS  - Extra definies needed by gromacs
#  GROMACS_PKG          - The name of the pkg-config package needed
#  GROMACS_VERSION      - Gromacs lib interface version
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

find_package(PkgConfig)
list(LENGTH GROMACS_FIND_COMPONENTS GROMACS_NUM_COMPONENTS_WANTED)
if(${GROMACS_NUM_COMPONENTS_WANTED} LESS 1)
  message(FATAL_ERROR "No gromacs component to search given")
elseif(${GROMACS_NUM_COMPONENTS_WANTED} GREATER 1)
  message(FATAL_ERROR "We only support finding one gromacs component at this point, go and implement it ;-)")
elseif(${GROMACS_FIND_COMPONENTS} MATCHES "^lib(gmx|gromacs)(_d)?$")
  if(NOT GROMACS_PKG_OVERWRITE)
    set(GROMACS_PKG "${GROMACS_FIND_COMPONENTS}")
  else()
    set(GROMACS_PKG "${GROMACS_PKG_OVERWRITE}")
  endif()
  string(REGEX REPLACE "^lib(.*)" "\\1" GROMACS_LIBRARY_NAME "${GROMACS_PKG}")
else()
  message(FATAL_ERROR "We do not support finding ${GROMACS_FIND_COMPONENTS}, go and implement it ;-)")
endif()

pkg_check_modules(PC_GROMACS ${GROMACS_PKG})

find_library(GROMACS_LIBRARY NAMES ${GROMACS_LIBRARY_NAME} HINTS ${PC_GROMACS_LIBRARY_DIRS} )
if (GROMACS_LIBRARY)
  include(CheckLibraryExists)
  check_library_exists("${GROMACS_LIBRARY}" GromacsVersion "" FOUND_GROMACS_VERSION)
  if(NOT FOUND_GROMACS_VERSION)
    message(FATAL_ERROR "Could not find GromacsVersion in ${GROMACS_LIBRARY}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set GROMACS_LIBRARY by hand which sets the gromacs lib and it's dependencies (i.e. -DGROMACS_LIBRARY='/path/to/libgmx.so;/path/to/libblas.so;/path/to/libm.so')!")
  endif(NOT FOUND_GROMACS_VERSION)
  check_library_exists("${GROMACS_LIBRARY}" init_mtop "" FOUND_GROMACS_INIT_MTOP)
  if(NOT FOUND_GROMACS_INIT_MTOP)
    message(FATAL_ERROR "Could not find init_mtop in the gromacs library, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. This most likely means that your gromacs version is too old, we need at least gromacs 4.0.7.") 
  endif(NOT FOUND_GROMACS_INIT_MTOP)
  set(GMX_VERSION 40)
  check_library_exists("${GROMACS_LIBRARY}" output_env_init "" FOUND_GROMACS_OUTPUT_ENV_INIT)
  if(FOUND_GROMACS_OUTPUT_ENV_INIT)
    set(GMX_VERSION 45)
  endif(FOUND_GROMACS_OUTPUT_ENV_INIT)
  check_library_exists("${GROMACS_LIBRARY}" init_domdec_vsites "" FOUND_GROMACS_INIT_DOMDEC_VSITES)
  if(FOUND_GROMACS_INIT_DOMDEC_VSITES)
    set(GMX_VERSION 50)
  endif(FOUND_GROMACS_INIT_DOMDEC_VSITES)
  set(GROMACS_VERSION ${GMX_VERSION})

  #Only set GROMACS_DEFINITIONS if GROMACS_LIBRARY was found
  if ("${GROMACS_LIBRARY}" MATCHES "lib[^/]*_d\\.[^.]*$")
    list(APPEND GMX_DEFS "-DGMX_DOUBLE")
  endif ("${GROMACS_LIBRARY}" MATCHES "lib[^/]*_d\\.[^.]*$")
  if (PC_GROMACS_CFLAGS_OTHER)
    foreach(DEF ${PC_GROMACS_CFLAGS_OTHER})
      if (${DEF} MATCHES "^-D")
        list(APPEND GMX_DEFS ${DEF})
      endif (${DEF} MATCHES "^-D")
    endforeach(DEF)
    list(REMOVE_DUPLICATES GMX_DEFS)
  endif (PC_GROMACS_CFLAGS_OTHER)
  set(GROMACS_DEFINITIONS "${GMX_DEFS}")

else(GROMACS_LIBRARY)
  set(GMX_VERSION 45)
endif (GROMACS_LIBRARY)

if ("${GROMACS_PKG}" MATCHES "libgmx")
  if (${GMX_VERSION} EQUAL 40)
    find_path(GROMACS_INCLUDE_DIR tpxio.h HINTS ${PC_GROMACS_INCLUDE_DIRS})
  else(${GMX_VERSION} EQUAL 40)
   find_path(GROMACS_INCLUDE_DIR gromacs/tpxio.h HINTS ${PC_GROMACS_INCLUDE_DIRS})
  endif(${GMX_VERSION} EQUAL 40)
elseif("${GROMACS_PKG}" MATCHES "libgromacs")
  find_path(GROMACS_INCLUDE_DIR gromacs/fileio/tpxio.h HINTS ${PC_GROMACS_INCLUDE_DIRS})
endif("${GROMACS_PKG}" MATCHES "libgmx")

set(GROMACS_LIBRARIES "${GROMACS_LIBRARY}" )
set(GROMACS_INCLUDE_DIRS ${GROMACS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GROMACS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GROMACS DEFAULT_MSG GROMACS_LIBRARY GROMACS_INCLUDE_DIR)

mark_as_advanced(GROMACS_INCLUDE_DIR GROMACS_LIBRARY GROMACS_PKG)
