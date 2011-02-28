# - Find libgmx
# Find the native libgmx headers and libraries.
#
#  GMX_INCLUDE_DIRS - where to find gromacs/tpxio.h, etc.
#  GMX_LIBRARIES    - List of libraries when using libgmx.
#  GMX_FOUND        - True if libgmx found.

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

pkg_check_modules(PC_GMX_D libgmx_d)
pkg_check_modules(PC_GMX libgmx)

#double over single
find_path(GMX_INCLUDE_DIR gromacs/tpxio.h HINTS ${PC_GMX_D_INCLUDE_DIRS} ${PC_GMX_INCLUDE_DIRS})
find_library(GMX_LIBRARY NAMES gmx_d gmx HINTS ${PC_GMX_D_LIBRARY_DIRS} ${PC_GMX_LIBRARY_DIRS})

set(GMX_LIBRARIES ${GMX_LIBRARY} )
set(GMX_INCLUDE_DIRS ${GMX_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GMX_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GMX DEFAULT_MSG GMX_LIBRARY GMX_INCLUDE_DIR )

mark_as_advanced(GMX_INCLUDE_DIR GMX_LIBRARY )
