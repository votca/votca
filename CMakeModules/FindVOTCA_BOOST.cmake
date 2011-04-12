# - Find libvotca_boost
# Find the native libvotca_boost headers and libraries.
#
#  VOTCA_BOOST_INCLUDE_DIRS - where to find headers etc.
#  VOTCA_BOOST_LIBRARIES    - List of libraries when using expat.
#  VOTCA_BOOST_FOUND        - True if expat found.
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

pkg_check_modules(PC_VOTCA_BOOST libvotca_boost)
find_path(VOTCA_BOOST_INCLUDE_DIR boost/algorithm/string/trim.hpp HINTS ${PC_VOTCA_BOOST_INCLUDE_DIRS} )
find_library(VOTCA_BOOST_LIBRARY NAMES votca_boost HINTS ${PC_VOTCA_BOOST_LIBRARY_DIRS} )

set(VOTCA_BOOST_LIBRARIES "${VOTCA_BOOST_LIBRARY}" )
set(VOTCA_BOOST_INCLUDE_DIRS "${VOTCA_BOOST_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_BOOST_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VOTCA_BOOST DEFAULT_MSG VOTCA_BOOST_LIBRARY VOTCA_BOOST_INCLUDE_DIR)

mark_as_advanced(VOTCA_BOOST_INCLUDE_DIR VOTCA_BOOST_LIBRARY )
