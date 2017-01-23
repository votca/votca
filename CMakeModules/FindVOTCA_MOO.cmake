# - Find libvotca_moo
# Find the native libvotca_moo headers and libraries.
#
#  VOTCA_MOO_INCLUDE_DIRS - where to find votca/moo/version.h, etc.
#  VOTCA_MOO_LIBRARIES    - List of libraries when using expat.
#  VOTCA_MOO_FOUND        - True if expat found.
#  VOTCA_MOO_HAS_SQLITE3  - True if votca moo was build with sqlite3 support
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

pkg_check_modules(PC_VOTCA_MOO libvotca_moo)
find_path(VOTCA_MOO_INCLUDE_DIR votca/moo/global.h HINTS ${PC_VOTCA_MOO_INCLUDE_DIRS})

find_library(VOTCA_MOO_LIBRARY NAMES votca_moo HINTS ${PC_VOTCA_MOO_LIBRARY_DIRS} )

if("${VOTCA_MOO_LIBRARY}" MATCHES "libvotca_moo[^;]*\\.a")
    if(PC_VOTCA_MOO_LIBRARIES)
      list(REMOVE_ITEM PC_VOTCA_MOO_LIBRARIES votca_moo)
      foreach (LIB ${PC_VOTCA_MOO_LIBRARIES})
        find_library(VOTCA_MOO_${LIB} NAMES ${LIB} HINTS ${PC_VOTCA_MOO_LIBRARY_DIRS} )
        list(APPEND VT_DEP_LIBRARIES ${VOTCA_MOO_${LIB}})
        unset(VOTCA_MOO_${LIB} CACHE)
      endforeach(LIB)
    endif(PC_VOTCA_MOO_LIBRARIES)
    set(VOTCA_MOO_DEP_LIBRARIES "${VT_DEP_LIBRARIES}" CACHE FILEPATH "votca moo depency libs (only needed for static (.a) libvotca_moo")
endif("${VOTCA_MOO_LIBRARY}" MATCHES "libvotca_moo[^;]*\\.a")

set(VOTCA_MOO_LIBRARIES "${VOTCA_MOO_LIBRARY};${VOTCA_MOO_DEP_LIBRARIES}" )
set(VOTCA_MOO_INCLUDE_DIRS "${VOTCA_MOO_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_MOO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VOTCA_MOO DEFAULT_MSG VOTCA_MOO_LIBRARY VOTCA_MOO_INCLUDE_DIR )

if (VOTCA_MOO_FOUND AND NOT VOTCA_MOO_LIBRARY STREQUAL "votca_moo")
  include(CheckLibraryExists)
  check_library_exists("${VOTCA_MOO_LIBRARY};${VOTCA_MOO_DEP_LIBRARIES}" VotcaMooFromC "" FOUND_VOTCA_MOO_VERSION)
  if(NOT FOUND_VOTCA_MOO_VERSION)
    message(FATAL_ERROR "Could not find VotcaMooFromC in ${VOTCA_MOO_LIBRARY};${VOTCA_MOO_DEP_LIBRARIES}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set VOTCA_MOO_LIBRARY and VOTCA_MOO_DEP_LIBRARIES by hand, which set votca_moo lib  it's depencies (i.e. -DVOTCA_MOO_LIBRARY='/path/to/libvotca_moo.so' -VOTCA_MOO_DEP_LIBRARIES='/path/to/libgsl.so;/path/to/libm.so') !")
  endif(NOT FOUND_VOTCA_MOO_VERSION)
endif ()

mark_as_advanced(VOTCA_MOO_INCLUDE_DIR VOTCA_MOO_LIBRARY )
