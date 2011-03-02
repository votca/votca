# - Find libgmx_d
# Find the native libgmx headers and libraries.
#
#  GMX_D_INCLUDE_DIRS - where to find gromacs/tpxio.h, etc.
#  GMX_D_LIBRARIES    - List of libraries when using libgmx_d.
#  GMX_D_FOUND        - True if libgmx_d found.
#  GMX_D_DEFINITIONS  - Extra definies needed by libgmx_d
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

pkg_check_modules(PC_GMX_D libgmx_d)
#avoid PC_GMX_D_LIBRARIES being empty
list(INSERT PC_GMX_D_LIBRARIES 0 gmx_d)
list(REMOVE_DUPLICATES PC_GMX_D_LIBRARIES)

find_path(GMX_D_INCLUDE_DIR gromacs/tpxio.h HINTS ${PC_GMX_D_INCLUDE_DIRS})
if (NOT GMX_D_LIBRARY)
  #add deps
  foreach (LIB ${PC_GMX_D_LIBRARIES})
    find_library(GMX_D_${LIB} NAMES ${LIB} HINTS ${PC_GMX_D_LIBRARY_DIRS} )
    list(APPEND GMX_D_LIBRARY ${GMX_D_${LIB}})
  endforeach(LIB)
endif (NOT GMX_D_LIBRARY)

if (NOT GMX_D_DEFINITIONS)
  list(INSERT GMX_D_DEFINITIONS 0 "-DGMX_DOUBLE")
  foreach(DEF ${PC_GMX_D_CFLAGS_OTHER})
    if (${DEF} MATCHES "^-D")
      list(APPEND GMX_D_DEFINITIONS ${DEF})
    endif (${DEF} MATCHES "^-D")
  endforeach(DEF)
  list(REMOVE_DUPLICATES GMX_D_DEFINITIONS)
endif (NOT GMX_D_DEFINITIONS)

set(GMX_D_LIBRARIES ${GMX_D_LIBRARY} )
set(GMX_D_NCLUDE_DIRS ${GMX_D_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GMX_D_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GMX_D DEFAULT_MSG GMX_D_LIBRARY GMX_D_INCLUDE_DIR )

mark_as_advanced(GMX_D_INCLUDE_DIR GMX_D_LIBRARY GMX_D_DEFINITIONS )
