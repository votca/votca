# - Find libgromacs
# Find the native libgromacs headers and libraries.
#
#  GROMACS_INCLUDE_DIRS - where to find gromacs/t, etc.
#  GROMACS_LIBRARIES    - List of libraries when using expat.
#  GROMACS_FOUND        - True if expat found.

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

pkg_check_modules(PC_GROMACS_D libgromacs_d)
pkg_check_modules(PC_GROMACS libgromacs)

#double over single
find_path(GROMACS_INCLUDE_DIR gromacs/legacyheaders/tpxio.h HINTS ${PC_GROMACS_D_INCLUDE_DIRS} ${PC_GROMACS_INCLUDE_DIRS})
find_library(GROMACS_LIBRARY NAMES gromacs_d gromacs HINTS ${PC_GROMACS_D_LIBRARY_DIRS} ${PC_GROMACS_LIBRARY_DIRS})

set(GROMACS_LIBRARIES ${GROMACS_LIBRARY} )
set(GROMACS_INCLUDE_DIRS ${GROMACS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GROMACS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GROMACS DEFAULT_MSG GROMACS_LIBRARY GROMACS_INCLUDE_DIR )

mark_as_advanced(GROMACS_INCLUDE_DIR GROMACS_LIBRARY )
