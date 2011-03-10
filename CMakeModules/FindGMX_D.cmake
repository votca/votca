# - Find libgmx_d
# Find the native libgmx_d headers and libraries.
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
find_path(GMX_D_INCLUDE_DIR gromacs/tpxio.h HINTS ${PC_GMX_D_INCLUDE_DIRS})
find_library(GMX_D_LIBRARY NAMES gmx_d HINTS ${PC_GMX_D_LIBRARY_DIRS} )

#add deps
if (NOT BUILD_SHARED_LIBS AND PC_GMX_D_LIBRARIES)
  list(REMOVE_ITEM PC_GMX_D_LIBRARIES gmx_d)
  foreach (LIB ${PC_GMX_D_LIBRARIES})
    find_library(GMX_D_${LIB} NAMES ${LIB} HINTS ${PC_GMX_D_LIBRARY_DIRS} )
    list(APPEND GMX_D_LIBRARY ${GMX_D_${LIB}})
    unset(GMX_D_${LIB} CACHE)
  endforeach(LIB)
  set(GMX_D_LIBRARY "${GMX_D_LIBRARY}" CACHE FILEPATH "GMX_D lib and it's deps" FORCE)
endif (NOT BUILD_SHARED_LIBS AND PC_GMX_D_LIBRARIES)

set(GMX_D_DEFINITIONS "" CACHE STRING "extra GMX_D definitions")
if (NOT GMX_D_DEFINITIONS)
  list(INSERT PC_GMX_D_CFLAGS_OTHER 0 "-DGMX_DOUBLE")
  foreach(DEF ${PC_GMX_D_CFLAGS_OTHER})
    if (${DEF} MATCHES "^-D")
      list(APPEND GMX_D_DEFINITIONS ${DEF})
    endif (${DEF} MATCHES "^-D")
  endforeach(DEF)
  list(REMOVE_DUPLICATES GMX_D_DEFINITIONS)
  set(GMX_D_DEFINITIONS "${GMX_D_DEFINITIONS}" CACHE STRING "extra GMX_D definitions" FORCE)
endif (NOT GMX_D_DEFINITIONS)

set(GMX_D_LIBRARIES ${GMX_D_LIBRARY} )
set(GMX_D_INCLUDE_DIRS ${GMX_D_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GMX_D_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GMX_D DEFAULT_MSG GMX_D_LIBRARY GMX_D_INCLUDE_DIR )

if (GMX_D_FOUND)
  include(CheckLibraryExists)
  check_library_exists("${GMX_D_LIBRARY}" GromacsVersion "" FOUND_GMX_D_VERSION)
  if(NOT FOUND_GMX_D_VERSION)
    message(FATAL_ERROR "Could not find GromacsVersion in ${GMX_D_LIBRARY}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set GMX_D_LIBRARY by hand and include all it's deps in there (i.e. -DGMX_D_LIBRARY='/path/to/libgmx_d.so;/path/to/libblas.so;/path/to/libm.so') !")
  endif(NOT FOUND_GMX_D_VERSION)
  check_library_exists("${GMX_D_LIBRARY}" init_mtop "" FOUND_GMX_D_INIT_MTOP)
  if(NOT FOUND_GMX_D_INIT_MTOP)
    message(FATAL_ERROR "Could not find GromacsVersion in ${GMX_D_LIBRARY}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. This mostly means that your libgmx_d is too old, we need at least libgmx_d-4.0.7.") 
  endif(NOT FOUND_GMX_D_INIT_MTOP)
endif (GMX_D_FOUND)

mark_as_advanced(GMX_D_INCLUDE_DIR GMX_D_LIBRARY GMX_D_DEFINITIONS )
