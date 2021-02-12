###############################################################################
# 
# CMake script for finding the Eigen library.
# 
# http://eigen.tuxfamily.org/index.php?title=Main_Page
# 
# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
# Redistribution and use is allowed according to the terms of the 2-clause BSD
# license.
# 
# 
# Input variables:
# 
# - Eigen_ROOT_DIR (optional): When specified, header files and libraries
#   will be searched for in `${Eigen_ROOT_DIR}/include` and
#   `${Eigen_ROOT_DIR}/libs` respectively, and the default CMake search order
#   will be ignored. When unspecified, the default CMake search order is used.
#   This variable can be specified either as a CMake or environment variable.
#   If both are set, preference is given to the CMake variable.
#   Use this variable for finding packages installed in a nonstandard location,
#   or for enforcing that one of multiple package installations is picked up.
# 
# Cache variables (not intended to be used in CMakeLists.txt files)
# 
# - Eigen_INCLUDE_DIR: Absolute path to package headers.
# 
# 
# Output variables:
# 
# - Eigen_FOUND: Boolean that indicates if the package was found
# - Eigen_INCLUDE_DIRS: Paths to the necessary header files
# - Eigen_VERSION: Version of Eigen library found
# - Eigen_DEFINITIONS: Definitions to be passed on behalf of eigen
# 
# 
# Example usage:
# 
#   # Passing the version means Eigen_FOUND will only be TRUE if a
#   # version >= the provided version is found.
#   find_package(Eigen 3.1.2)
#   if(NOT Eigen_FOUND)
#     # Error handling
#   endif()
#   ...
#   add_definitions(${Eigen_DEFINITIONS})
#   ...
#   include_directories(${Eigen_INCLUDE_DIRS} ...)
# 
###############################################################################

find_package(PkgConfig)
pkg_check_modules(PC_EIGEN eigen3)


find_package(Eigen3 QUIET NO_MODULE PATHS ${PC_EIGEN_PREFIX}/ NO_DEFAULT_PATH)
if (NOT EIGEN3_FOUND)
  # Search user-installed locations first, so that we prefer user installs
  # to system installs where both exist.
  list(APPEND EIGEN3_CHECK_INCLUDE_DIRS
    /usr/local/include
    /usr/local/homebrew/include # Mac OS X
    /opt/local/var/macports/software # Mac OS X.
    /opt/local/include
    /usr/include)
  # Additional suffixes to try appending to each search path.
  list(APPEND EIGEN3_CHECK_PATH_SUFFIXES
    eigen3 # Default root directory for Eigen.
    Eigen/include/eigen3 # Windows (for C:/Program Files prefix) < 3.3
    Eigen3/include/eigen3 ) # Windows (for C:/Program Files prefix) >= 3.3

  # Search supplied hint directories first if supplied.
  find_path(EIGEN3_INCLUDE_DIR
    NAMES Eigen/Core
    HINTS ${EIGEN3_INCLUDE_DIR_HINTS}
    PATHS ${EIGEN3_CHECK_INCLUDE_DIRS}
    PATH_SUFFIXES ${EIGEN3_CHECK_PATH_SUFFIXES})

  if (NOT EIGEN3_INCLUDE_DIR OR
      NOT EXISTS ${EIGEN3_INCLUDE_DIR})
    message(ERROR 
      "Could not find eigen3 include directory, set EIGEN_INCLUDE_DIR to "
      "path to eigen3 include directory, e.g. /usr/local/include/eigen3.")
  endif (NOT EIGEN3_INCLUDE_DIR OR
    NOT EXISTS ${EIGEN3_INCLUDE_DIR})

  # Mark internally as found, then verify. EIGEN_REPORT_NOT_FOUND() unsets
  # if called.
  set(EIGEN3_FOUND TRUE)
endif()



if (EIGEN3_INCLUDE_DIR)
  set(EIGEN3_VERSION_FILE ${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h)
  if (NOT EXISTS ${EIGEN3_VERSION_FILE})
    message(ERROR
      "Could not find file: ${EIGEN3_VERSION_FILE} "
      "containing version information in Eigen install located at: "
      "${EIGEN3_INCLUDE_DIR}.")
  else (NOT EXISTS ${EIGEN3_VERSION_FILE})
    file(READ ${EIGEN3_VERSION_FILE} EIGEN3_VERSION_FILE_CONTENTS)

    string(REGEX MATCH "#define EIGEN_WORLD_VERSION [0-9]+"
      EIGEN3_WORLD_VERSION "${EIGEN3_VERSION_FILE_CONTENTS}")
    string(REGEX REPLACE "#define EIGEN_WORLD_VERSION ([0-9]+)" "\\1"
      EIGEN3_WORLD_VERSION "${EIGEN3_WORLD_VERSION}")

    string(REGEX MATCH "#define EIGEN_MAJOR_VERSION [0-9]+"
      EIGEN3_MAJOR_VERSION "${EIGEN3_VERSION_FILE_CONTENTS}")
    string(REGEX REPLACE "#define EIGEN_MAJOR_VERSION ([0-9]+)" "\\1"
      EIGEN3_MAJOR_VERSION "${EIGEN3_MAJOR_VERSION}")

    string(REGEX MATCH "#define EIGEN_MINOR_VERSION [0-9]+"
      EIGEN3_MINOR_VERSION "${EIGEN3_VERSION_FILE_CONTENTS}")
    string(REGEX REPLACE "#define EIGEN_MINOR_VERSION ([0-9]+)" "\\1"
      EIGEN3_MINOR_VERSION "${EIGEN3_MINOR_VERSION}")

    # This is on a single line s/t CMake does not interpret it as a list of
    # elements and insert ';' separators which would result in 3.;2.;0 nonsense.
    set(EIGEN3_VERSION ${EIGEN3_WORLD_VERSION}.${EIGEN3_MAJOR_VERSION}.${EIGEN3_MINOR_VERSION})
  endif (NOT EXISTS ${EIGEN3_VERSION_FILE})
endif (EIGEN3_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EIGEN3 REQUIRED_VARS EIGEN3_INCLUDE_DIR VERSION_VAR EIGEN3_VERSION)

mark_as_advanced(EIGEN3_INCLUDE_DIR EIGEN3_VERSION)



