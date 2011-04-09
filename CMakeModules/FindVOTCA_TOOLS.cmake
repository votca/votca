# - Find libgmx
# Find the native libgmx headers and libraries.
#
#  VOTCA_TOOLS_INCLUDE_DIRS - where to find votca/tools/version.h, etc.
#  VOTCA_TOOLS_LIBRARIES    - List of libraries when using expat.
#  VOTCA_TOOLS_FOUND        - True if expat found.
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

pkg_check_modules(PC_VOTCA_TOOLS libvotca_tools)
find_path(VOTCA_TOOLS_INCLUDE_DIR votca/tools/version.h HINTS ${PC_VOTCA_TOOLS_INCLUDE_DIRS})
find_path(VOTCA_BOOST_INCLUDE_DIR boost/algorithm/string/trim.hpp HINTS ${PC_VOTCA_TOOLS_INCLUDE_DIRS} ${VOTCA_TOOLS_INCLUDE_DIR}/votca )
find_library(VOTCA_TOOLS_LIBRARY NAMES votca_tools HINTS ${PC_VOTCA_TOOLS_LIBRARY_DIRS} )

#add deps
if (NOT BUILD_SHARED_LIBS AND PC_VOTCA_TOOLS_LIBRARIES)
  list(REMOVE_ITEM PC_VOTCA_TOOLS_LIBRARIES votca_tools)
  foreach (LIB ${PC_VOTCA_TOOLS_LIBRARIES})
    find_library(VOTCA_TOOLS_${LIB} NAMES ${LIB} HINTS ${PC_VOTCA_TOOLS_LIBRARY_DIRS} )
    list(APPEND VOTCA_TOOLS_LIBRARY ${VOTCA_TOOLS_${LIB}})
    unset(VOTCA_TOOLS_${LIB} CACHE)
  endforeach(LIB)
  set(VOTCA_TOOLS_LIBRARY "${VOTCA_TOOLS_LIBRARY}" CACHE FILEPATH "VOTCA_TOOLS lib and it's deps" FORCE)
endif (NOT BUILD_SHARED_LIBS AND PC_VOTCA_TOOLS_LIBRARIES)

set(VOTCA_TOOLS_LIBRARIES "${VOTCA_TOOLS_LIBRARY}" )
set(VOTCA_TOOLS_INCLUDE_DIRS "${VOTCA_TOOLS_INCLUDE_DIR};${VOTCA_BOOST_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_TOOLS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VOTCA_TOOLS DEFAULT_MSG VOTCA_TOOLS_LIBRARY VOTCA_TOOLS_INCLUDE_DIR VOTCA_BOOST_INCLUDE_DIR)

if (VOTCA_TOOLS_FOUND)
  include(CheckVotcaToolsExists)
  check_votca_tools_exists("${VOTCA_TOOLS_LIBRARY}" ToolsVersionStr "" FOUND_VOTCA_TOOLS_VERSION)
  if(NOT FOUND_VOTCA_TOOLS_VERSION)
    message(FATAL_ERROR "Could not find votca::tools::ToolsVersionStr in ${VOTCA_TOOLS_LIBRARY}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set VOTCA_TOOLS_LIBRARY by hand and include all it's deps in there (i.e. -DVOTCA_TOOLS_LIBRARY='/path/to/libvotca_tools.so;/path/to/libgsl.so;/path/to/libm.so') !")
  endif(NOT FOUND_VOTCA_TOOLS_VERSION)
endif (VOTCA_TOOLS_FOUND)

mark_as_advanced(VOTCA_TOOLS_INCLUDE_DIR VOTCA_TOOLS_LIBRARY )
