# - Find libgromacs_d
# Find the native libgromacs_d headers and libraries.
#
#  GROMACS_D_INCLUDE_DIRS - where to find gromacs/tpxio.h, etc.
#  GROMACS_D_LIBRARIES    - List of libraries when using libgromacs_d.
#  GROMACS_D_FOUND        - True if libgromacs_d found.
#  GROMACS_D_DEFINITIONS  - Extra definies needed by libgromacs_d
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

pkg_check_modules(PC_GROMACS_D libgromacs_d)
find_path(GROMACS_D_INCLUDE_DIR gromacs/legacyheaders/tpxio.h HINTS ${PC_GROMACS_D_INCLUDE_DIRS})
find_library(GROMACS_D_LIBRARY NAMES gromacs_d HINTS ${PC_GROMACS_D_LIBRARY_DIRS} )

#add deps
if (NOT BUILD_SHARED_LIBS AND PC_GROMACS_D_LIBRARIES)
  list(REMOVE_ITEM PC_GROMACS_D_LIBRARIES gromacs_d)
  foreach (LIB ${PC_GROMACS_D_LIBRARIES})
    find_library(GROMACS_D_${LIB} NAMES ${LIB} HINTS ${PC_GROMACS_D_LIBRARY_DIRS} )
    list(APPEND GROMACS_D_LIBRARY ${GROMACS_D_${LIB}})
    unset(GROMACS_D_${LIB} CACHE)
  endforeach(LIB)
  set(GROMACS_D_LIBRARY "${GROMACS_D_LIBRARY}" CACHE FILEPATH "GROMACS_D lib and it's deps" FORCE)
endif (NOT BUILD_SHARED_LIBS AND PC_GROMACS_D_LIBRARIES)

set(GROMACS_D_DEFINITIONS "" CACHE STRING "extra GROMACS_D definitions")
if (NOT GROMACS_D_DEFINITIONS)
  list(INSERT PC_GROMACS_D_CFLAGS_OTHER 0 "-DGMX_DOUBLE")
  foreach(DEF ${PC_GROMACS_D_CFLAGS_OTHER})
    if (${DEF} MATCHES "^-D")
      list(APPEND GROMACS_D_DEFINITIONS ${DEF})
    endif (${DEF} MATCHES "^-D")
  endforeach(DEF)
  list(REMOVE_DUPLICATES GROMACS_D_DEFINITIONS)
  set(GROMACS_D_DEFINITIONS "${GROMACS_D_DEFINITIONS}" CACHE STRING "extra GROMACS_D definitions" FORCE)
endif (NOT GROMACS_D_DEFINITIONS)

set(GROMACS_D_LIBRARIES ${GROMACS_D_LIBRARY} )
set(GROMACS_D_INCLUDE_DIRS ${GROMACS_D_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GROMACS_D_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GROMACS_D DEFAULT_MSG GROMACS_D_LIBRARY GROMACS_D_INCLUDE_DIR )

if (GROMACS_D_FOUND)
  include(CheckLibraryExists)
  check_library_exists("${GROMACS_D_LIBRARY}" GromacsVersion "" FOUND_GROMACS_D_VERSION)
  if(NOT FOUND_GROMACS_D_VERSION)
    message(FATAL_ERROR "Could not find GromacsVersion in ${GROMACS_D_LIBRARY}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set GROMACS_D_LIBRARY by hand and include all it's deps in there (i.e. -DGROMACS_D_LIBRARY='/path/to/libgromacs_d.so;/path/to/libblas.so;/path/to/libm.so') !")
  endif(NOT FOUND_GROMACS_D_VERSION)
endif (GROMACS_D_FOUND)

mark_as_advanced(GROMACS_D_INCLUDE_DIR GROMACS_D_LIBRARY GROMACS_D_DEFINITIONS )
