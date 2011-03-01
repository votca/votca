# - Find libgromacs
# Find the native libgromacs headers and libraries.
#
#  GROMACS_INCLUDE_DIRS - where to find gromacs/tpxio.h, etc.
#  GROMACS_LIBRARIES    - List of libraries when using libgromacs.
#  GROMACS_FOUND        - True if libgromacs found.
#  GROMACS_PKG          - Whether we are using libgromacs oder libgromacs_d
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

pkg_check_modules(PC_GROMACS libgromacs_d)
set(GROMACS_PKG "libgromacs_d")
find_path(GROMACS_INCLUDE_DIR gromacs/legacyheaders/tpxio.h HINTS ${PC_GROMACS_INCLUDE_DIRS})
foreach (LIB ${PC_GROMACS_LIBRARIES})
  find_library(GROMACS_${LIB} NAMES ${LIB} HINTS ${PC_GROMACS_LIBRARY_DIRS} )
  list(APPEND GROMACS_LIBRARY ${GROMACS_${LIB}})
endforeach(LIB)
if ("${GROMACS_LIBRARY}" STREQUAL "NOTFOUND" OR "${GROMACS_INCLUDE_DIR}" STREQUAL "NOTFOUND")
  pkg_check_modules(PC_GROMACS libgromacs)
  find_path(GROMACS_INCLUDE_DIR gromacs/legacyheaders/tpxio.h HINTS ${PC_GROMACS_INCLUDE_DIRS})
  foreach (LIB ${PC_GROMACS_LIBRARIES})
    find_library(GROMACS_${LIB} NAMES ${LIB} HINTS ${PC_GROMACS_LIBRARY_DIRS} )
    list(APPEND GROMACS_LIBRARY ${GROMACS_${LIB}})
  endforeach(LIB)
  set(GROMACS_PKG "libgromacs")
endif ("${GROMACS_LIBRARY}" STREQUAL "NOTFOUND" OR "${GROMACS_INCLUDE_DIR}" STREQUAL "NOTFOUND")

foreach(DEF ${PC_GROMACS_CFLAGS_OTHER})
  if (${DEF} MATCHES "^-D")
    list(APPEND GROMACS_DEFINITIONS ${DEF})
  endif (${DEF} MATCHES "^-D")
endforeach(DEF)

set(GROMACS_LIBRARIES ${GROMACS_LIBRARY} )
set(GROMACS_INCLUDE_DIRS ${GROMACS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GROMACS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GROMACS DEFAULT_MSG GROMACS_LIBRARY GROMACS_INCLUDE_DIR )

mark_as_advanced(GROMACS_INCLUDE_DIR GROMACS_LIBRARY GROMACS_DEFINITIONS )
