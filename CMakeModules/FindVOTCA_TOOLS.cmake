# - Find libgmx
# Find the native libgmx headers and libraries.
#
#  VOTCA_TOOLS_INCLUDE_DIRS - where to find votca/tools/version.h, etc.
#  VOTCA_TOOLS_LIBRARIES    - List of libraries when using expat.
#  VOTCA_TOOLS_FOUND        - True if expat found.

#=============================================================================
# Copyright 2011 VOTCA Development Team
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

find_package(PkgConfig)

pkg_check_modules(PC_VOTCA_TOOLS libvotca_tools)

#double over single
find_path(VOTCA_TOOLS_INCLUDE_DIR votca/tools/version.h HINTS ${PC_VOTCA_TOOLS_INCLUDE_DIRS})
find_path(VOTCA_BOOST_INCLUDE_DIR boost/algorithm/string/trim.hpp HINTS ${PC_VOTCA_TOOLS_INCLUDE_DIRS})
find_library(VOTCA_TOOLS_LIBRARY NAMES votca_tools HINTS ${PC_VOTCA_TOOLS_LIBRARY_DIRS})

set(VOTCA_TOOLS_LIBRARIES ${VOTCA_TOOLS_LIBRARY} )
set(VOTCA_TOOLS_INCLUDE_DIRS "${VOTCA_TOOLS_INCLUDE_DIR};${VOTCA_BOOST_INCLUDE_DIR}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VOTCA_TOOLS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VOTCA_TOOLS DEFAULT_MSG VOTCA_TOOLS_LIBRARY VOTCA_TOOLS_INCLUDE_DIR )
find_package_handle_standard_args(VOTCA_TOOLS DEFAULT_MSG VOTCA_TOOLS_LIBRARY VOTCA_BOOST_INCLUDE_DIR )

mark_as_advanced(VOTCA_TOOLS_INCLUDE_DIR VOTCA_BOOST_INCLUDE_DIR VOTCA_TOOLS_LIBRARY )
