# - Find libsqlite3
# Find the native libsqlite3 headers and libraries.
#
#  SQLITE3_INCLUDE_DIRS - where to find sqlite3.h, etc
#  SQLITE3_LIBRARIES    - List of libraries when using sqlite3.
#  SQLITE3_FOUND        - True if sqlite3 found.
#
# Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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

pkg_check_modules(PC_SQLITE3 sqlite3)

find_path(SQLITE3_INCLUDE_DIR sqlite3.h HINTS ${PC_SQLITE3_INCLUDE_DIRS})
find_library(SQLITE3_LIBRARY NAMES sqlite3 HINTS ${PC_SQLITE3_LIBRARY_DIRS} )

set(SQLITE3_LIBRARIES "${SQLITE3_LIBRARY}" )
set(SQLITE3_INCLUDE_DIRS "${SQLITE3_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SQLITE3 DEFAULT_MSG SQLITE3_LIBRARY SQLITE3_INCLUDE_DIR )

mark_as_advanced(SQLITE3_INCLUDE_DIR SQLITE3_LIBRARY )
