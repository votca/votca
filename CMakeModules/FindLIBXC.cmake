# - Find libvotca_moo
# Find the native libvotca_moo headers and libraries.
#
#  VOTCA_MOO_INCLUDE_DIRS - where to find votca/moo/version.h, etc.
#  VOTCA_MOO_LIBRARIES    - List of libraries when using expat.
#  VOTCA_MOO_FOUND        - True if expat found.
#  VOTCA_MOO_HAS_SQLITE3  - True if votca moo was build with sqlite3 support
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

# check if info is available via PkgConfig
find_package(PkgConfig)
pkg_check_modules(PC_libxc libxc)

find_path(LIBXC_INCLUDE_DIR NAMES xc.h HINTS HINTS ${PC_LIBXC_INCLUDE_DIRS})
find_library(LIBXC_LIBRARY NAMES xc HINTS ${PC_LIBXC_LIBRARY_DIRS} )


set(LIBXC_LIBRARIES "${LIBXC_LIBRARY}" )
set(LIBXC_INCLUDE_DIRS "${LIBXC_INCLUDE_DIR}" )

#include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_MOO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LIBXC DEFAULT_MSG LIBXC_LIBRARY LIBXC_INCLUDE_DIR )


mark_as_advanced(LIBXC_INCLUDE_DIR LIBXC_LIBRARY )
