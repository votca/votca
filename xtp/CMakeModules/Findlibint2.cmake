# - Find libint2
# Find libint2 through pkg-config
#
# Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
#
#  LIBINT2_INCLUDE_DIRS - where to find libint2 headers.
#  LIBINT2_LIBRARIES    - List of libraries when used by libint2.
#  LIBINT2_FOUND        - True if all libint2 componets were found.
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

function(_LIBINT2_GET_VERSION _OUT_ver _version_hdr)
  if(NOT EXISTS ${_version_hdr})
    message(FATAL_ERROR "Header file ${_version_hdr} not found (check value of LIBINT2_INCLUDE_DIR)")
  endif()
  file(STRINGS ${_version_hdr} _contents REGEX "#define LIBINT_VERSION[ \t]+")
  message(STATUS ${_contents})
  if(_contents)
    string(REGEX REPLACE ".*#define LIBINT_VERSION[ \t]+\"([0-9.]+)\".*" "\\1" ${_OUT_ver} "${_contents}")
    if(NOT ${${_OUT_ver}} MATCHES "^[0-9.]+$")
      message(FATAL_ERROR "Version parsing failed for LIBINT_VERSION in ${_version_hdr}, excepted a number but got ${${_OUT_ver}}!")
    endif()
    set(${_OUT_ver} ${${_OUT_ver}} PARENT_SCOPE)
  else()
    message(FATAL_ERROR "No LIBINT_VERSION  line found in include file ${_version_hdr}")
  endif()
endfunction()

find_package(PkgConfig)
if(COMMAND set_package_properties)
  set_package_properties(PkgConfig PROPERTIES TYPE RECOMMENDED PURPOSE "Used to detect libint2 package")
endif()

pkg_check_modules(PC_LIBINT2 libint2)
find_path(LIBINT2_INCLUDE_DIR libint2.hpp HINTS ${PC_LIBINT2_INCLUDE_DIRS})
if(LIBINT2_INCLUDE_DIR)
  _LIBINT2_GET_VERSION(LIBINT2_VERSION ${LIBINT2_INCLUDE_DIR}/libint2/config.h)
else()
  set(LIBINT2_VERSION 0)
endif()

find_library(LIBINT2_LIBRARY NAMES int2 HINTS ${PC_LIBINT2_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(libint2 REQUIRED_VARS LIBINT2_INCLUDE_DIR LIBINT2_LIBRARY VERSION_VAR LIBINT2_VERSION) 

if(libint2_FOUND)
  set(LIBINT2_LIBRARIES ${LIBINT2_LIBRARY})
  set(LIBINT2_INCLUDE_DIRS ${LIBINT2_INCLUDE_DIR})

  if(NOT Libint2::int2)
    add_library(Libint2::int2 UNKNOWN IMPORTED)
    set_target_properties(Libint2::int2 PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "CPP"
      IMPORTED_LOCATION "${LIBINT2_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${LIBINT2_INCLUDE_DIRS}")
  endif()
endif()
