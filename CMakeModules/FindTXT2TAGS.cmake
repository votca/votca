# Copyright (C) 2011 Votca Development Team
#
# This file was derived from FindGnuplot.cmake shipped with CMake 2.6.3.
#
# - this module looks for txt2tags
#
# Once done this will define
#
#  TXT2TAGS_FOUND - system has txt2tags
#  TXT2TAGS_EXECUTABLE - the txt2tags executable
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

FIND_PROGRAM(TXT2TAGS_EXECUTABLE
  NAMES 
  txt2tags
  txt2tags-2.5
  txt2tags-2.6
  PATHS
  ${CYGWIN_INSTALL_PATH}/bin
)

# handle the QUIETLY and REQUIRED arguments and set TXT2TAGS_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(TXT2TAGS DEFAULT_MSG TXT2TAGS_EXECUTABLE)

IF(NOT TXT2TAGS_FOUND)
  message("txt2tags not found, help cmake to find it by setting TXT2TAGS_EXECUTABLE")
ENDIF(NOT TXT2TAGS_FOUND)

MARK_AS_ADVANCED( TXT2TAGS_EXECUTABLE )

