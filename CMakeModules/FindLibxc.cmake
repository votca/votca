# - Find libxc
# Find libxc through pkg-config
#
# Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
#
#  Libxc_INCLUDE_DIRS - where to find libxc headers.
#  Libxc_LIBRARIES    - List of libraries when used by libxc.
#  Libxc_FOUND        - True if all libxc componets were found.
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
if(COMMAND set_package_properties)
  set_package_properties(PkgConfig PROPERTIES TYPE RECOMMENDED PURPOSE "Used to detect libxc package")
endif()

pkg_check_modules(PC_Libxc libxc)
find_path(Libxc_INCLUDE_DIR xc.h HINTS ${PC_Libxc_INCLUDE_DIRS})

find_library(Libxc_LIBRARY NAMES xc HINTS ${PC_Libxc_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(Libxc REQUIRED_VARS Libxc_INCLUDE_DIR Libxc_LIBRARY) 

if(Libxc_FOUND)
  set(Libxc_LIBRARIES ${Libxc_LIBRARY})
  set(Libxc_INCLUDE_DIRS ${Libxc_INCLUDE_DIR})

  if(NOT Libxc::xc)
    add_library(Libxc::xc UNKNOWN IMPORTED)
    set_target_properties(Libxc::xc PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      IMPORTED_LOCATION "${Libxc_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${Libxc_INCLUDE_DIRS}")
  endif()
endif()
