# - Find gsl
# Find the native GSL headers and libraries.
#
#  GSL_INCLUDE_DIRS - where to find gsl/gsl_linalg.h, etc.
#  GSL_LIBRARIES    - List of libraries when using gsl.
#  GSL_FOUND        - True if gsl found.
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

pkg_check_modules(PC_GSL gsl)
find_path(GSL_INCLUDE_DIR gsl/gsl_linalg.h HINTS ${PC_GSL_INCLUDE_DIRS})
find_library(GSL_LIBRARY NAMES gsl HINTS ${PC_GSL_LIBRARY_DIRS} )

#added gsl depends
if (NOT BUILD_SHARED_LIBS)
  list(REMOVE_ITEM PC_GSL_LIBRARIES gsl)
  foreach (LIB ${PC_GSL_LIBRARIES})
    find_library(GSL_${LIB} NAMES ${LIB} HINTS ${PC_GSL_LIBRARY_DIRS} )
    list(APPEND GSL_LIBRARY ${GSL_${LIB}})
    unset(GSL_${LIB} CACHE)
  endforeach(LIB)
endif (NOT BUILD_SHARED_LIBS)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GSL DEFAULT_MSG GSL_LIBRARY GSL_INCLUDE_DIR )

set(GSL_LIBRARIES ${GSL_LIBRARY} )
set(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR} )

mark_as_advanced(GSL_INCLUDE_DIR GSL_LIBRARY)
