# - Find libgromacs_d
# Find the native libgromacs headers and libraries.
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
#avoid PC_GROMACS_D_LIBRARIES being empty
list(INSERT PC_GROMACS_D_LIBRARIES 0 gromacs_d)
list(REMOVE_DUPLICATES PC_GROMACS_D_LIBRARIES)

find_path(GROMACS_D_INCLUDE_DIR gromacs/legacyheaders/tpxio.h HINTS ${PC_GROMACS_D_INCLUDE_DIRS})
if (NOT GROMACS_D_LIBRARY)
  #add deps
  foreach (LIB ${PC_GROMACS_D_LIBRARIES})
    find_library(GROMACS_D_${LIB} NAMES ${LIB} HINTS ${PC_GROMACS_D_LIBRARY_DIRS} )
    list(APPEND GROMACS_D_LIBRARY ${GROMACS_D_${LIB}})
  endforeach(LIB)
endif (NOT GROMACS_D_LIBRARY)

if (NOT GROMACS_D_DEFINITIONS)
  list(INSERT GROMACS_D_DEFINITIONS 0 "-DGMX_DOUBLE")
  foreach(DEF ${PC_GROMACS_D_CFLAGS_OTHER})
    if (${DEF} MATCHES "^-D")
      list(APPEND GROMACS_D_DEFINITIONS ${DEF})
    endif (${DEF} MATCHES "^-D")
  endforeach(DEF)
  list(REMOVE_DUPLICATES GROMACS_D_DEFINITIONS)
endif (NOT GROMACS_D_DEFINITIONS)

set(GROMACS_D_LIBRARIES ${GROMACS_D_LIBRARY} )
set(GROMACS_D_INCLUDE_DIRS ${GROMACS_D_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GROMACS_D_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GROMACS_D DEFAULT_MSG GROMACS_D_LIBRARY GROMACS_D_INCLUDE_DIR )

mark_as_advanced(GROMACS_D_INCLUDE_DIR GROMACS_D_LIBRARY GROMACS_D_DEFINITIONS )
