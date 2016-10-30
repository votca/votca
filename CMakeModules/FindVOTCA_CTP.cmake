# - Find libvotca_ctp
# Find the native libvotca_ctp headers and libraries.
#
#  VOTCA_CTP_INCLUDE_DIRS - where to find votca/ctp/version.h, etc.
#  VOTCA_CTP_LIBRARIES    - List of libraries when using expat.
#  VOTCA_CTP_FOUND        - True if expat found.
#  VOTCA_CTP_HAS_SQLITE3  - True if votca ctp was build with sqlite3 support
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

pkg_check_modules(PC_VOTCA_CTP libvotca_ctp)
find_path(VOTCA_CTP_INCLUDE_DIR votca/ctp/version.h HINTS ${PC_VOTCA_CTP_INCLUDE_DIRS})

find_library(VOTCA_CTP_LIBRARY NAMES votca_ctp HINTS ${PC_VOTCA_CTP_LIBRARY_DIRS} )

if("${VOTCA_CTP_LIBRARY}" MATCHES "libvotca_ctp[^;]*\\.a")
    if(PC_VOTCA_CTP_LIBRARIES)
      list(REMOVE_ITEM PC_VOTCA_CTP_LIBRARIES votca_ctp)
      foreach (LIB ${PC_VOTCA_CTP_LIBRARIES})
        find_library(VOTCA_CTP_${LIB} NAMES ${LIB} HINTS ${PC_VOTCA_CTP_LIBRARY_DIRS} )
        list(APPEND VT_DEP_LIBRARIES ${VOTCA_CTP_${LIB}})
        unset(VOTCA_CTP_${LIB} CACHE)
      endforeach(LIB)
    endif(PC_VOTCA_CTP_LIBRARIES)
    set(VOTCA_CTP_DEP_LIBRARIES "${VT_DEP_LIBRARIES}" CACHE FILEPATH "votca ctp depency libs (only needed for static (.a) libvotca_ctp")
endif("${VOTCA_CTP_LIBRARY}" MATCHES "libvotca_ctp[^;]*\\.a")

set(VOTCA_CTP_LIBRARIES "${VOTCA_CTP_LIBRARY};${VOTCA_CTP_DEP_LIBRARIES}" )
set(VOTCA_CTP_INCLUDE_DIRS "${VOTCA_CTP_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_CTP_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VOTCA_CTP DEFAULT_MSG VOTCA_CTP_LIBRARY VOTCA_CTP_INCLUDE_DIR )

if (VOTCA_CTP_FOUND AND NOT VOTCA_CTP_LIBRARY STREQUAL "votca_ctp")
  include(CheckLibraryExists)
  check_library_exists("${VOTCA_CTP_LIBRARY};${VOTCA_CTP_DEP_LIBRARIES}" VotcaMd2QmFromC "" FOUND_VOTCA_CTP_VERSION)
  if(NOT FOUND_VOTCA_CTP_VERSION)
    message(FATAL_ERROR "Could not find VotcaMd2QmFromC in ${VOTCA_CTP_LIBRARY};${VOTCA_CTP_DEP_LIBRARIES}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set VOTCA_CTP_LIBRARY and VOTCA_CTP_DEP_LIBRARIES by hand, which set votca_ctp lib  it's depencies (i.e. -DVOTCA_CTP_LIBRARY='/path/to/libvotca_ctp.so' -VOTCA_CTP_DEP_LIBRARIES='/path/to/libgsl.so;/path/to/libm.so') !")
  endif(NOT FOUND_VOTCA_CTP_VERSION)
endif ()

mark_as_advanced(VOTCA_CTP_INCLUDE_DIR VOTCA_CTP_LIBRARY )
