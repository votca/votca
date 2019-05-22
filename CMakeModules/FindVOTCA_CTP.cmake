# - Find libvotca_ctp
# Find the native libvotca_ctp headers and libraries.
#
#  VOTCA_CTP_INCLUDE_DIRS - where to find votca/ctp/version.h, etc.
#  VOTCA_CTP_LIBRARIES    - List of libraries when using expat.
#  VOTCA_CTP_FOUND        - True if expat found.
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

pkg_check_modules(PC_VOTCA_CTP libvotca_ctp)
find_path(VOTCA_CTP_INCLUDE_DIR votca/ctp/ctpapplication.h HINTS ${PC_VOTCA_CTP_INCLUDE_DIRS})

find_library(VOTCA_CTP_LIBRARY NAMES votca_ctp HINTS ${PC_VOTCA_CTP_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_CTP_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VOTCA_CTP DEFAULT_MSG VOTCA_CTP_LIBRARY VOTCA_CTP_INCLUDE_DIR )

if (VOTCA_CTP_FOUND AND NOT VOTCA_CTP_LIBRARY STREQUAL "votca_ctp")
  include(CheckLibraryExists)
  check_library_exists("${VOTCA_CTP_LIBRARY};${VOTCA_CTP_DEP_LIBRARIES}" VotcaMd2QmFromC "" FOUND_VOTCA_CTP_VERSION)
  if(NOT FOUND_VOTCA_CTP_VERSION)
    message(FATAL_ERROR "Could not find VotcaMd2QmFromC in ${VOTCA_CTP_LIBRARY};${VOTCA_CTP_DEP_LIBRARIES}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set VOTCA_CTP_LIBRARY and VOTCA_CTP_DEP_LIBRARIES by hand, which set votca_ctp lib  it's depencies (i.e. -DVOTCA_CTP_LIBRARY='/path/to/libvotca_ctp.so' -VOTCA_CTP_DEP_LIBRARIES='/path/to/libgsl.so;/path/to/libm.so') !")
  endif(NOT FOUND_VOTCA_CTP_VERSION)
endif ()

if(VOTCA_CTP_FOUND)
  set(VOTCA_CTP_LIBRARIES "${VOTCA_CTP_LIBRARY}" )
  set(VOTCA_CTP_INCLUDE_DIRS "${VOTCA_CTP_INCLUDE_DIR}" )

  if(NOT TARGET VOTCA::votca_ctp)
    add_library(VOTCA::votca_ctp UNKNOWN IMPORTED)
    set_target_properties(VOTCA::votca_ctp PROPERTIES
      IMPORTED_LOCATION "${VOTCA_CTP_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${VOTCA_CTP_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(VOTCA_CTP_INCLUDE_DIR VOTCA_CTP_LIBRARY )
