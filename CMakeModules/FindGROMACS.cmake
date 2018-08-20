# - Finds parts of gromacs
# Find the native gromacs compents headers and libraries.
#
#  GROMACS_INCLUDE_DIRS - where to find gromacs headers.
#  GROMACS_LIBRARIES    - List of libraries when used by gromacs.
#  GROMACS_FOUND        - True if all gromacs componets were found.
#  GROMACS_DEFINITIONS  - Extra definies needed by gromacs
#  GROMACS_VERSION      - Gromacs lib interface version
#
# Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

function(_GROMACS_GET_VERSION _OUT_ver _version_hdr)
  file(STRINGS ${_version_hdr} _contents REGEX "#define GMX_VERSION[ \t]+")
  if(_contents)
    string(REGEX REPLACE ".*#define GMX_VERSION[ \t]+([0-9.]+).*" "\\1" ${_OUT_ver} "${_contents}")
    if(("${${_OUT_ver}}" STREUQAL "") OR (NOT ${${_OUT_ver}} MATCHES "[0-9]+"))
        message(FATAL_ERROR "Version parsing failed for GMX_VERSION in ${_version_hdr}!")
    endif()
    set(${_OUT_ver} ${${_OUT_ver}} PARENT_SCOPE)
 elseif(EXISTS ${_version_hdr})
    message(FATAL_ERROR "No GMX_VERSION in ${_version_hdr}")
 else()
     message(FATAL_ERROR "No GMX_VERSION line found in include file ${_version_hdr}")
  endif()
endfunction()

find_package(PkgConfig)
pkg_check_modules(PC_GROMACS_D libgromacs_d)
pkg_check_modules(PC_GROMACS libgromacs)

find_library(GROMACS_LIBRARY NAMES gromacs_d gromacs HINTS ${PC_GROMACS_D_LIBRARY_DIRS} ${PC_GROMACS_LIBRARY_DIRS})
if (GROMACS_LIBRARY AND NOT GROMACS_LIBRARY STREQUAL "gromacs")
  include(CheckCXXLibraryExists)
  check_cxx_library_exists("${GROMACS_LIBRARY}" gmx_version "" FOUND_GROMACS_VERSION_CXX)
  if(NOT FOUND_GROMACS_VERSION_CXX)
    message(FATAL_ERROR "Could not find a suitable gromacs library. gmx_version is not  defined in the gromacs library, that is very very strange, take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. This most likely means that your gromacs version is too old, we need at least gromacs 2016!")
  endif()
  check_cxx_library_exists("${GROMACS_LIBRARY}" gmx_is_single_precision "" FOUND_GMX_IS_SINGLE_PRECISION)
  check_cxx_library_exists("${GROMACS_LIBRARY}" gmx_is_double_precision "" FOUND_GMX_IS_DOUBLE_PRECISION)
  if(FOUND_GMX_IS_DOUBLE_PRECISION)
    set(GROMACS_DEFINITIONS "-DGMX_DOUBLE=1")
  elseif(FOUND_GMX_IS_SINGLE_PRECISION)
    set(GROMACS_DEFINITIONS "-DGMX_DOUBLE=0")
  elseif(NOT FOUND_GMX_IS_SINGLE_PRECISION)
    message(FATAL_ERROR "Could not find a suitable gromacs library. Neither gmx_is_single_precision nor gmx_is_double_precision is defined in the gromacs library, that is very very strange, take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. This most likely means that your gromacs version is too old, we need at least gromacs 2016!")
  endif()
endif()

find_path(GROMACS_INCLUDE_DIR gromacs/version.h HINTS ${PC_GROMACS_D_INCLUDE_DIRS} ${PC_GROMACS_INCLUDE_DIRS})
if(GROMACS_VERSION)
elseif(GROMACS_INCLUDE_DIR AND EXISTS ${GROMACS_INCLUDE_DIR}/gromacs/version.h)
  _GROMACS_GET_VERSION(GROMACS_VERSION ${GROMACS_INCLUDE_DIR}/gromacs/version.h)
else()
  set(GROMACS_VERSION 0)
endif()

set(GROMACS_LIBRARIES "${GROMACS_LIBRARY}" )
set(GROMACS_INCLUDE_DIRS "${GROMACS_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GROMACS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GROMACS REQUIRED_VARS GROMACS_LIBRARY GROMACS_INCLUDE_DIR VERSION_VAR GROMACS_VERSION)

mark_as_advanced(GROMACS_LIBRARY GROMACS_INCLUDE_DIR GROMACS_VERSION)
