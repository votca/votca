# - Find libvotca_kmcold
# Find the native libvotca_kmcold headers and libraries.
#
#  VOTCA_KMCOLD_INCLUDE_DIRS - where to find votca/kmcold/version.h, etc.
#  VOTCA_KMCOLD_LIBRARIES    - List of libraries when using expat.
#  VOTCA_KMCOLD_FOUND        - True if expat found.
#  VOTCA_KMCOLD_HAS_SQLITE3  - True if votca kmcold was build with sqlite3 support
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

pkg_check_modules(PC_VOTCA_KMCOLD libvotca_kmcold)
find_path(VOTCA_KMCOLD_INCLUDE_DIR votca/kmcold/kmc.h HINTS ${PC_VOTCA_KMCOLD_INCLUDE_DIRS})

find_library(VOTCA_KMCOLD_LIBRARY NAMES votca_kmcold HINTS ${PC_VOTCA_KMCOLD_LIBRARY_DIRS} )

if("${VOTCA_KMCOLD_LIBRARY}" MATCHES "libvotca_kmcold[^;]*\\.a")
    if(PC_VOTCA_KMCOLD_LIBRARIES)
      list(REMOVE_ITEM PC_VOTCA_KMCOLD_LIBRARIES votca_kmcold)
      foreach (LIB ${PC_VOTCA_KMCOLD_LIBRARIES})
        find_library(VOTCA_KMCOLD_${LIB} NAMES ${LIB} HINTS ${PC_VOTCA_KMCOLD_LIBRARY_DIRS} )
        list(APPEND VT_DEP_LIBRARIES ${VOTCA_KMCOLD_${LIB}})
        unset(VOTCA_KMCOLD_${LIB} CACHE)
      endforeach(LIB)
    endif(PC_VOTCA_KMCOLD_LIBRARIES)
    set(VOTCA_KMCOLD_DEP_LIBRARIES "${VT_DEP_LIBRARIES}" CACHE FILEPATH "votca kmcold depency libs (only needed for static (.a) libvotca_kmcold")
endif("${VOTCA_KMCOLD_LIBRARY}" MATCHES "libvotca_kmcold[^;]*\\.a")

set(VOTCA_KMCOLD_LIBRARIES "${VOTCA_KMCOLD_LIBRARY};${VOTCA_KMCOLD_DEP_LIBRARIES}" )
set(VOTCA_KMCOLD_INCLUDE_DIRS "${VOTCA_KMCOLD_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_KMCOLD_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VOTCA_KMCOLD DEFAULT_MSG VOTCA_KMCOLD_LIBRARY VOTCA_KMCOLD_INCLUDE_DIR )

if (VOTCA_KMCOLD_FOUND)
  include(CheckLibraryExists)
  check_library_exists("${VOTCA_KMCOLD_LIBRARY};${VOTCA_KMCOLD_DEP_LIBRARIES}" VotcaKmcOldFromC "" FOUND_VOTCA_KMCOLD_VERSION)
  if(NOT FOUND_VOTCA_KMCOLD_VERSION)
    message(FATAL_ERROR "Could not find VotcaKmcOldFromC in ${VOTCA_KMCOLD_LIBRARY};${VOTCA_KMCOLD_DEP_LIBRARIES}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set VOTCA_KMCOLD_LIBRARY and VOTCA_KMCOLD_DEP_LIBRARIES by hand, which set votca_kmcold lib  it's depencies (i.e. -DVOTCA_KMCOLD_LIBRARY='/path/to/libvotca_kmcold.so" -VOTCA_KMCOLD_DEP_LIBRARIES="/path/to/libgsl.so;/path/to/libm.so') !")
  endif(NOT FOUND_VOTCA_KMCOLD_VERSION)
endif (VOTCA_KMCOLD_FOUND)

mark_as_advanced(VOTCA_KMCOLD_INCLUDE_DIR VOTCA_KMCOLD_LIBRARY )
