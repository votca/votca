# - Find libgromacs
# Find the native libgromacs headers and libraries.
#
#  GROMACS_INCLUDE_DIRS - where to find gromacs/tpxio.h, etc.
#  GROMACS_LIBRARIES    - List of libraries when using libgromacs.
#  GROMACS_FOUND        - True if libgromacs found.
#  GROMACS_DEFINITIONS  - Extra definies needed by libgromacs
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

pkg_check_modules(PC_GROMACS libgromacs)
find_path(GROMACS_INCLUDE_DIR gromacs/legacyheaders/tpxio.h HINTS ${PC_GROMACS_INCLUDE_DIRS})
find_library(GROMACS_LIBRARY NAMES gromacs HINTS ${PC_GROMACS_LIBRARY_DIRS} )

#add deps
if (NOT BUILD_SHARED_LIBS AND PC_GROMACS_LIBRARIES)
  list(REMOVE_ITEM PC_GROMACS_LIBRARIES gromacs)
  foreach (LIB ${PC_GROMACS_LIBRARIES})
    find_library(GROMACS_${LIB} NAMES ${LIB} HINTS ${PC_GROMACS_LIBRARY_DIRS} )
    list(APPEND GROMACS_LIBRARY ${GROMACS_${LIB}})
    unset(GROMACS_${LIB} CACHE)
  endforeach(LIB)
  set(GROMACS_LIBRARY "${GROMACS_LIBRARY}" CACHE FILEPATH "GROMACS lib and it's deps" FORCE)
endif (NOT BUILD_SHARED_LIBS AND PC_GROMACS_LIBRARIES)

set(GROMACS_DEFINITIONS "" CACHE STRING "extra GROMACS definitions")
if (NOT GROMACS_DEFINITIONS AND PC_GROMACS_CFLAGS_OTHER)
  foreach(DEF ${PC_GROMACS_CFLAGS_OTHER})
    if (${DEF} MATCHES "^-D")
      list(APPEND GROMACS_DEFINITIONS ${DEF})
    endif (${DEF} MATCHES "^-D")
  endforeach(DEF)
  set(GROMACS_DEFINITIONS "${GROMACS_DEFINITIONS}" CACHE STRING "extra GROMACS definitions" FORCE)
endif (NOT GROMACS_DEFINITIONS AND PC_GROMACS_CFLAGS_OTHER)

set(GROMACS_LIBRARIES ${GROMACS_LIBRARY} )
set(GROMACS_INCLUDE_DIRS ${GROMACS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GROMACS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GROMACS DEFAULT_MSG GROMACS_LIBRARY GROMACS_INCLUDE_DIR )

if (GROMACS_FOUND)
  include(CheckLibraryExists)
  check_library_exists("${GROMACS_LIBRARY}" GromacsVersion "" FOUND_GROMACS_VERSION)
  if(NOT FOUND_GROMACS_VERSION)
    message(FATAL_ERROR "Could not find GromacsVersion in ${GROMACS_LIBRARY}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set GROMACS_LIBRARY by hand and include all it's deps in there (i.e. -DGROMACS_LIBRARY='/path/to/libgromacs.so;/path/to/libblas.so;/path/to/libm.so') !")
  endif(NOT FOUND_GROMACS_VERSION)
endif (GROMACS_FOUND)

mark_as_advanced(GROMACS_INCLUDE_DIR GROMACS_LIBRARY GROMACS_DEFINITIONS )
