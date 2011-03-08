# - Find libgmx
# Find the native libgmx headers and libraries.
#
#  GMX_INCLUDE_DIRS - where to find gromacs/tpxio.h, etc.
#  GMX_LIBRARIES    - List of libraries when using libgmx.
#  GMX_FOUND        - True if libgmx found.
#  GMX_DEFINITIONS  - Extra definies needed by libgmx
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

pkg_check_modules(PC_GMX libgmx)
find_path(GMX_INCLUDE_DIR gromacs/tpxio.h HINTS ${PC_GMX_INCLUDE_DIRS})
find_library(GMX_LIBRARY NAMES gmx HINTS ${PC_GMX_LIBRARY_DIRS} )

#add deps
if (NOT BUILD_SHARED_LIBS)
  list(REMOVE_ITEM PC_GMX_LIBRARIES gmx)
  foreach (LIB ${PC_GMX_LIBRARIES})
    find_library(GMX_${LIB} NAMES ${LIB} HINTS ${PC_GMX_LIBRARY_DIRS} )
    list(APPEND GMX_LIBRARY ${GMX_${LIB}})
    unset(GMX_${LIB} CACHE)
  endforeach(LIB)
endif (NOT BUILD_SHARED_LIBS)

set(GMX_DEFINITIONS "" CACHE STRING "extra GMX DEFINITIONS")
if (NOT GMX_DEFINITIONS)
  foreach(DEF ${PC_GMX_CFLAGS_OTHER})
    if (${DEF} MATCHES "^-D")
      list(APPEND GMX_DEFINITIONS ${DEF})
    endif (${DEF} MATCHES "^-D")
  endforeach(DEF)
endif (NOT GMX_DEFINITIONS)

set(GMX_LIBRARIES ${GMX_LIBRARY} )
set(GMX_INCLUDE_DIRS ${GMX_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GMX_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GMX DEFAULT_MSG GMX_LIBRARY GMX_INCLUDE_DIR )

if (GMX_FOUND)
  include(CheckLibraryExists)
  check_library_exists("${GMX_LIBRARY}" GromacsVersion "" FOUND_GMX_VERSION)
  if(NOT FOUND_GMX_VERSION)
    message(FATAL_ERROR "Could not find GromacsVersion in ${GMX_LIBRARY}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set GMX_LIBRARY by hand and include all it's deps in there (i.e. -DGMX_LIBRARY='/path/to/libgmx.so;/path/to/libblas.so;/path/to/libm.so') !")
  endif(NOT FOUND_GMX_VERSION)
  check_library_exists("${GMX_LIBRARY}" init_mtop "" FOUND_GMX_INIT_MTOP)
  if(NOT FOUND_GMX_INIT_MTOP)
    message(FATAL_ERROR "Could not find GromacsVersion in ${GMX_LIBRARY}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. This mostly means that your libgmx version is too old, we need at least libgmx-4.0.7.") 
  endif(NOT FOUND_GMX_INIT_MTOP)
endif (GMX_FOUND)

mark_as_advanced(GMX_INCLUDE_DIR GMX_LIBRARY GMX_DEFINITIONS )
