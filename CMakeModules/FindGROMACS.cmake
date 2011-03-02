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
#avoid PC_GROMACS_LIBRARIES being empty
list(INSERT PC_GROMACS_LIBRARIES 0 gromacs)
list(REMOVE_DUPLICATES PC_GROMACS_LIBRARIES)

find_path(GROMACS_INCLUDE_DIR gromacs/legacyheaders/tpxio.h HINTS ${PC_GROMACS_INCLUDE_DIRS})
if (NOT GROMACS_LIBRARY)
  #add deps
  foreach (LIB ${PC_GROMACS_LIBRARIES})
    find_library(GROMACS_${LIB} NAMES ${LIB} HINTS ${PC_GROMACS_LIBRARY_DIRS} )
    list(APPEND GROMACS_LIBRARY ${GROMACS_${LIB}})
  endforeach(LIB)
endif (NOT GROMACS_LIBRARY)

if (NOT GROMACS_DEFINITIONS)
  foreach(DEF ${PC_GROMACS_CFLAGS_OTHER})
      if (${DEF} MATCHES "^-D")
        list(APPEND GROMACS_DEFINITIONS ${DEF})
      endif (${DEF} MATCHES "^-D")
    endforeach(DEF)
endif (NOT GROMACS_DEFINITIONS)

set(GROMACS_LIBRARIES ${GROMACS_LIBRARY} )
set(GROMACS_NCLUDE_DIRS ${GROMACS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GROMACS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GROMACS DEFAULT_MSG GROMACS_LIBRARY GROMACSD_INCLUDE_DIR )

mark_as_advanced(GROMACS_INCLUDE_DIR GROMACS_LIBRARY GROMACS_DEFINITIONS )
