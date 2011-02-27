# - Find gsl
# Find the native GSL headers and libraries.
#
#  GSL_INCLUDE_DIRS - where to find expat.h, etc.
#  GSL_LIBRARIES    - List of libraries when using expat.
#  GSL_FOUND        - True if expat found.

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

pkg_check_modules(PC_GSL gsl)
find_path(GSL_INCLUDE_DIR gsl/gsl_linalg.h HINTS ${PC_GSL_INCLUDE_DIRS})

find_library(GSL_LIBRARY NAMES gsl HINTS ${PC_GSL_LIBRARY_DIRS} )
find_library(GSL_CBLAS_LIBRARY NAMES gslcblas HINTS ${PC_GSL_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GSL DEFAULT_MSG GSL_LIBRARY GSL_INCLUDE_DIR )
find_package_handle_standard_args(GSL DEFAULT_MSG GSL_CBLAS_LIBRARY GSL_INCLUDE_DIR )

set(GSL_LIBRARIES "${GSL_LIBRARY};${GSL_CBLAS_LIBRARY}" )
set(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR} )

mark_as_advanced(GSL_INCLUDE_DIR GSL_LIBRARY GSL_CBLAS_LIBRARY)
