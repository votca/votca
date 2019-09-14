# Copyright (C) 2019 Votca Development Team
#
# This file was derived from FindGnuplot.cmake shipped with CMake 2.6.3.
#
# - this module looks for inkscape
#
# Once done this will define
#
#  INKSCAPE_FOUND - system has inkscape
#  INKSCAPE_EXECUTABLE - the inkscape executable
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
#=============================================================================
# Copyright 2002-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

INCLUDE(FindCygwin)

FIND_PROGRAM(INKSCAPE_EXECUTABLE
  NAMES
  inkscape
  PATHS
  ${CYGWIN_INSTALL_PATH}/bin
)

if(INKSCAPE_EXECUTABLE)
  execute_process(
    COMMAND "${INKSCAPE_EXECUTABLE}" --version
    OUTPUT_VARIABLE _INKSCAPE_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _inkscape_version_result
  )
  if(_inkscape_version_result)
    message(WARNING "Unable to determine inkscape version: ${_inkscape_version_result}")
  endif()
  string(REGEX REPLACE "^Inkscape ([0-9.]+) \\(.*$" "\\1" INKSCAPE_VERSION "${_INKSCAPE_VERSION}")
  if(NOT INKSCAPE_VERSION MATCHES "[0-9.]+")
    message(WARNING "Unable to parse inkscape version: ${_INKSCAPE_VERSION}")
  endif()
  if(INKSCAPE_VERSION VERSION_LESS "1.0")
    set(INKSCAPE_EXPORT_FLAG "-E")
  else()
    set(INKSCAPE_EXPORT_FLAG "--export-file")
  endif()
else()
  message(WARNING "inkscape not found, help cmake to find it by setting INKSCAPE_EXECUTABLE")
endif()

# handle the QUIETLY and REQUIRED arguments and set INKSCAPE_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(INKSCAPE REQUIRED_VARS INKSCAPE_EXECUTABLE INKSCAPE_EXPORT_FLAG VERSION_VAR INKSCAPE_VERSION)

mark_as_advanced( INKSCAPE_EXECUTABLE )

