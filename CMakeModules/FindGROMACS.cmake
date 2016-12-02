# - Finds parts of gromacs
# Find the native gromacs compents headers and libraries.
#
#  GROMACS_INCLUDE_DIRS - where to find gromacs headers.
#  GROMACS_LIBRARIES    - List of libraries when used by gromacs.
#  GROMACS_FOUND        - True if all gromacs componets were found.
#  GROMACS_DEFINITIONS  - Extra definies needed by gromacs
#  GROMACS_VERSION      - Gromacs lib interface version
#
# Copyright 2009-2015 The VOTCA Development Team (http://www.votca.org)
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
pkg_check_modules(PC_GROMACS_D libgromacs_d)
pkg_check_modules(PC_GROMACS libgromacs)

find_library(GROMACS_LIBRARY NAMES gromacs_d gromacs HINTS ${PC_GROMACS_D_LIBRARY_DIRS} ${PC_GROMACS_LIBRARY_DIRS})
if (GROMACS_LIBRARY)
  include(CheckLibraryExists)
  include(CheckCXXLibraryExists)
  check_library_exists("${GROMACS_LIBRARY}" GromacsVersion "" FOUND_GROMACS_VERSION)
  if(NOT FOUND_GROMACS_VERSION)
    check_cxx_library_exists("${GROMACS_LIBRARY}" gmx_version "" FOUND_GROMACS_VERSION_CXX)
  endif()
  if(NOT FOUND_GROMACS_VERSION AND NOT FOUND_GROMACS_VERSION_CXX)
    message(FATAL_ERROR "Could not find GromacsVersion in ${GROMACS_LIBRARY}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set GROMACS_LIBRARY by hand which sets the gromacs lib and it's dependencies (i.e. -DGROMACS_LIBRARY='/path/to/libgmx.so;/path/to/libblas.so;/path/to/libm.so')!")
  endif()
  check_library_exists("${GROMACS_LIBRARY}" gmx_gpu_sharing_supported "" FOUND_GROMACS_GMX_GPU_SHARING_SUPPORTED)
  #check is above
  if(FOUND_GROMACS_VERSION_CXX)
    set(GROMACS_VERSION 52)
  elseif(FOUND_GROMACS_GMX_GPU_SHARING_SUPPORTED)
    set(GROMACS_VERSION 51)
  else()
    message(FATAL_ERROR "Could not find gmx_version nor gmx_gpu_sharing_supported in the gromacs library, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. This most likely means that your gromacs version is too old, we need at least gromacs 5.1 !") 
  endif()
  check_cxx_library_exists("${GROMACS_LIBRARY}" gmx_is_single_precision "" FOUND_GMX_IS_SINGLE_PRECISION)
  check_cxx_library_exists("${GROMACS_LIBRARY}" gmx_is_double_precision "" FOUND_GMX_IS_DOUBLE_PRECISION)
  if(FOUND_GMX_IS_DOUBLE_PRECISION AND GROMACS_VERSION GREATER 51)
    set(GROMACS_DEFINITIONS "-DGMX_DOUBLE=1")
  elseif(FOUND_GMX_IS_SINGLE_PRECISION AND GROMACS_VERSION GREATER 51)
    set(GROMACS_DEFINITIONS "-DGMX_DOUBLE=0")
  elseif(FOUND_GMX_IS_DOUBLE_PRECISION)
    set(GROMACS_DEFINITIONS "-DGMX_DOUBLE")
  elseif(NOT FOUND_GMX_IS_SINGLE_PRECISION)
    message(FATAL_ERROR "Could not find neither gmx_is_single_precision nor gmx_is_double_precision in the gromacs library, that is very very strange, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. This most likely means that your gromacs version is too old, we need at least gromacs 5.1 !") 
  endif()
endif (GROMACS_LIBRARY)

find_path(GROMACS_INCLUDE_DIR gromacs/fileio/tpxio.h HINTS ${PC_GROMACS_D_INCLUDE_DIRS} ${PC_GROMACS_INCLUDE_DIRS})

set(GROMACS_LIBRARIES "${GROMACS_LIBRARY}" )
set(GROMACS_INCLUDE_DIRS "${GROMACS_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GROMACS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GROMACS DEFAULT_MSG GROMACS_LIBRARY GROMACS_INCLUDE_DIR GROMACS_VERSION)

mark_as_advanced(GROMACS_LIBRARY GROMACS_INCLUDE_DIR GROMACS_VERSION)
