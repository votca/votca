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
find_path(VOTCA_BOOST_INCLUDE_DIR boost/algorithm/string/trim.hpp HINTS ${PC_VOTCA_TOOLS_INCLUDE_DIRS})
list(APPEND VOTCA_TOOLS_INCLUDE_DIR ${VOTCA_BOOST_INCLUDE_DIR})

list(INSERT PC_VOTCA_TOOLS_LIBRARIES 0 votca_tools)
list(REMOVE_DUPLICATES PC_VOTCA_TOOLS_LIBRARIES)
foreach (LIB ${PC_VOTCA_TOOLS_LIBRARIES})
  find_library(VOTCA_TOOLS_${LIB} NAMES ${LIB} HINTS ${PC_VOTCA_TOOLS_LIBRARY_DIRS} )
  list(APPEND VOTCA_TOOLS_LIBRARY ${VOTCA_TOOLS_${LIB}})
endforeach(LIB)

set(VOTCA_TOOLS_LIBRARIES "${VOTCA_TOOLS_LIBRARY}" )
set(VOTCA_TOOLS_INCLUDE_DIRS "${VOTCA_TOOLS_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_TOOLS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VOTCA_TOOLS DEFAULT_MSG VOTCA_TOOLS_LIBRARY VOTCA_TOOLS_INCLUDE_DIR )

mark_as_advanced(VOTCA_TOOLS_INCLUDE_DIR VOTCA_TOOLS_LIBRARY )
