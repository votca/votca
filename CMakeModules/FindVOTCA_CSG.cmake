# - Find libgmx
# Find the native libgmx headers and libraries.
#
#  VOTCA_CSG_INCLUDE_DIRS - where to find votca/csg/version.h, etc.
#  VOTCA_CSG_LIBRARIES    - List of libraries when using expat.
#  VOTCA_CSG_FOUND        - True if expat found.
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

pkg_check_modules(PC_VOTCA_CSG libvotca_csg)
find_path(VOTCA_CSG_INCLUDE_DIR votca/csg/version.h HINTS ${PC_VOTCA_CSG_INCLUDE_DIRS})
list(APPEND VOTCA_CSG_INCLUDE_DIR ${VOTCA_BOOST_INCLUDE_DIR})

list(INSERT PC_VOTCA_CSG_LIBRARIES 0 votca_csg)
list(REMOVE_DUPLICATES PC_VOTCA_CSG_LIBRARIES)
foreach (LIB ${PC_VOTCA_CSG_LIBRARIES})
  find_library(VOTCA_CSG_${LIB} NAMES ${LIB} HINTS ${PC_VOTCA_CSG_LIBRARY_DIRS} )
  list(APPEND VOTCA_CSG_LIBRARY ${VOTCA_CSG_${LIB}})
endforeach(LIB)

set(VOTCA_CSG_LIBRARIES "${VOTCA_CSG_LIBRARY}" )
set(VOTCA_CSG_INCLUDE_DIRS "${VOTCA_CSG_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_CSG_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VOTCA_CSG DEFAULT_MSG VOTCA_CSG_LIBRARY VOTCA_CSG_INCLUDE_DIR )

mark_as_advanced(VOTCA_CSG_INCLUDE_DIR VOTCA_CSG_LIBRARY )
