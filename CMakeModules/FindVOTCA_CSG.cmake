# - Find libvotca_csg
# Find the native libvotca_csg headers and libraries.
#
#  VOTCA_CSG_INCLUDE_DIRS - where to find votca/csg/version.h, etc.
#  VOTCA_CSG_LIBRARIES    - List of libraries when using expat.
#  VOTCA_CSG_FOUND        - True if expat found.
#  VOTCA_CSG_HAS_SQLITE3  - True if votca csg was build with sqlite3 support
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

pkg_check_modules(PC_VOTCA_CSG libvotca_csg)
find_path(VOTCA_CSG_INCLUDE_DIR votca/csg/version.h HINTS ${PC_VOTCA_CSG_INCLUDE_DIRS})

find_library(VOTCA_CSG_LIBRARY NAMES votca_csg HINTS ${PC_VOTCA_CSG_LIBRARY_DIRS} )

if("${VOTCA_CSG_LIBRARY}" MATCHES "libvotca_csg[^;]*\\.a")
    if(PC_VOTCA_CSG_LIBRARIES)
      list(REMOVE_ITEM PC_VOTCA_CSG_LIBRARIES votca_csg)
      foreach (LIB ${PC_VOTCA_CSG_LIBRARIES})
        find_library(VOTCA_CSG_${LIB} NAMES ${LIB} HINTS ${PC_VOTCA_CSG_LIBRARY_DIRS} )
        list(APPEND VT_DEP_LIBRARIES ${VOTCA_CSG_${LIB}})
        unset(VOTCA_CSG_${LIB} CACHE)
      endforeach(LIB)
    endif(PC_VOTCA_CSG_LIBRARIES)
    set(VOTCA_CSG_DEP_LIBRARIES "${VT_DEP_LIBRARIES}" CACHE FILEPATH "votca csg depency libs (only needed for static (.a) libvotca_csg")
endif("${VOTCA_CSG_LIBRARY}" MATCHES "libvotca_csg[^;]*\\.a")

set(VOTCA_CSG_LIBRARIES "${VOTCA_CSG_LIBRARY};${VOTCA_CSG_DEP_LIBRARIES}" )
set(VOTCA_CSG_INCLUDE_DIRS "${VOTCA_CSG_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_CSG_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VOTCA_CSG DEFAULT_MSG VOTCA_CSG_LIBRARY VOTCA_CSG_INCLUDE_DIR )

if (VOTCA_CSG_FOUND)
  include(CheckLibraryExists)
  check_library_exists("${VOTCA_CSG_LIBRARY};${VOTCA_CSG_DEP_LIBRARIES}" VotcaCsgFromC "" FOUND_VOTCA_CSG_VERSION)
  if(NOT FOUND_VOTCA_CSG_VERSION)
    message(FATAL_ERROR "Could not find VotcaCsgFromC in ${VOTCA_CSG_LIBRARY};${VOTCA_CSG_DEP_LIBRARIES}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set VOTCA_CSG_LIBRARY and VOTCA_CSG_DEP_LIBRARIES by hand, which set votca_csg lib  it's depencies (i.e. -DVOTCA_CSG_LIBRARY='/path/to/libvotca_csg.so" -VOTCA_CSG_DEP_LIBRARIES="/path/to/libgsl.so;/path/to/libm.so') !")
  endif(NOT FOUND_VOTCA_CSG_VERSION)
endif (VOTCA_CSG_FOUND)

mark_as_advanced(VOTCA_CSG_INCLUDE_DIR VOTCA_CSG_LIBRARY )
