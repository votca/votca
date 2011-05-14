# - Find libgmx
# Find the native libgmx headers and libraries.
#
#  VOTCA_KMCOLD_INCLUDE_DIRS - where to find votca/kmc/version.h, etc.
#  VOTCA_KMCOLD_LIBRARIES    - List of libraries when using expat.
#  VOTCA_KMCOLD_FOUND        - True if expat found.
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

pkg_check_modules(PC_VOTCA_KMCOLD libvotca_kmcold)
find_path(VOTCA_KMCOLD_INCLUDE_DIR votca/kmcold/kmc.h HINTS ${PC_VOTCA_KMCOLD_INCLUDE_DIRS})
list(APPEND VOTCA_KMCOLD_INCLUDE_DIR ${VOTCA_BOOST_INCLUDE_DIR})

list(INSERT PC_VOTCA_KMCOLD_LIBRARIES 0 votca_kmcold)
list(REMOVE_DUPLICATES PC_VOTCA_KMCOLD_LIBRARIES)
foreach (LIB ${PC_VOTCA_KMCOLD_LIBRARIES})
  find_library(VOTCA_KMCOLD_${LIB} NAMES ${LIB} HINTS ${PC_VOTCA_KMCOLD_LIBRARY_DIRS} )
  list(APPEND VOTCA_KMCOLD_LIBRARY ${VOTCA_KMCOLD_${LIB}})
endforeach(LIB)

set(VOTCA_KMCOLD_LIBRARIES "${VOTCA_KMCOLD_LIBRARY}" )
set(VOTCA_KMCOLD_INCLUDE_DIRS "${VOTCA_KMCOLD_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_KMCOLD_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VOTCA_KMCOLD DEFAULT_MSG VOTCA_KMCOLD_LIBRARY VOTCA_KMCOLD_INCLUDE_DIR )

mark_as_advanced(VOTCA_KMCOLD_INCLUDE_DIR VOTCA_KMCOLD_LIBRARY )
