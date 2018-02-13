# - Find libint
# Find the native libint headers and libraries.
#
#  LIBINT_INCLUDE_DIRS - where to find libint2.hpp, etc.
#  LIBINT_LIBRARIES    - List of libraries when using expat.
#  LIBINT_FOUND        - True if expat found.
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
pkg_check_modules(PC_libint2 libint2)

find_path(LIBINT_INCLUDE_DIR NAMES libint2.hpp HINTS HINTS ${PC_libint2_INCLUDE_DIRS})
find_library(LIBINT_LIBRARY NAMES int2 HINTS ${PC_libint2_LIBRARY_DIRS} )


set(LIBINT_LIBRARIES "${LIBINT_LIBRARY}" )
set(LIBINT_INCLUDE_DIRS "${LIBINT_INCLUDE_DIR}" )

#include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBINT_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LIBINT DEFAULT_MSG LIBINT_LIBRARY LIBINT_INCLUDE_DIR )


mark_as_advanced(LIBINT_INCLUDE_DIR LIBINT_LIBRARY )
