# Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
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

# try to find Apple Accelarate framework
find_library(ACCELERATE Accelerate REQUIRED)

# Eigen Overload needs also standalone LAPACKE
if (ACCELERATE)

    # Try to find LAPACKE lib
    #######################

    # Looking for include
    # -------------------

    # Add system include paths to search include
    # ------------------------------------------
    unset(_inc_env)
    set(ENV_LAPACKE_DIR "$ENV{LAPACKE_DIR}")
    set(ENV_LAPACKE_INCDIR "$ENV{LAPACKE_INCDIR}")
    if(ENV_LAPACKE_INCDIR)
      list(APPEND _inc_env "${ENV_LAPACKE_INCDIR}")
    elseif(ENV_LAPACKE_DIR)
      list(APPEND _inc_env "${ENV_LAPACKE_DIR}")
      list(APPEND _inc_env "${ENV_LAPACKE_DIR}/include")
      list(APPEND _inc_env "${ENV_LAPACKE_DIR}/include/lapacke")
    else()
      if(WIN32)
        string(REPLACE ":" ";" _inc_env "$ENV{INCLUDE}")
      else()
        string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
        list(APPEND _inc_env "${_path_env}")
      endif()
    endif()
    list(APPEND _inc_env "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")
    list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
    list(REMOVE_DUPLICATES _inc_env)


    # Try to find the lapacke header in the given paths
    # -------------------------------------------------
    # call cmake macro to find the header path
    if(LAPACKE_INCDIR)
      set(LAPACKE_lapacke.h_DIRS "LAPACKE_lapacke.h_DIRS-NOTFOUND")
      find_path(LAPACKE_lapacke.h_DIRS
        NAMES lapacke.h
        HINTS ${LAPACKE_INCDIR})
    else()
      if(LAPACKE_DIR)
        set(LAPACKE_lapacke.h_DIRS "LAPACKE_lapacke.h_DIRS-NOTFOUND")
        find_path(LAPACKE_lapacke.h_DIRS
          NAMES lapacke.h
          HINTS ${LAPACKE_DIR}
          PATH_SUFFIXES "include" "include/lapacke")
      else()
        set(LAPACKE_lapacke.h_DIRS "LAPACKE_lapacke.h_DIRS-NOTFOUND")
        find_path(LAPACKE_lapacke.h_DIRS
          NAMES lapacke.h
          HINTS ${_inc_env})
      endif()
    endif()
    mark_as_advanced(LAPACKE_lapacke.h_DIRS)

    # If found, add path to cmake variable
    # ------------------------------------
    if (LAPACKE_lapacke.h_DIRS)
      set(LAPACKE_INCLUDE_DIRS "${LAPACKE_lapacke.h_DIRS}")
    else ()
      set(LAPACKE_INCLUDE_DIRS "LAPACKE_INCLUDE_DIRS-NOTFOUND")
      if(NOT LAPACKE_FIND_QUIETLY)
        message(STATUS "Looking for lapacke -- lapacke.h not found")
      endif()
    endif()


    # Looking for lib
    # ---------------

    # Add system library paths to search lib
    # --------------------------------------
    unset(_lib_env)
    set(ENV_LAPACKE_LIBDIR "$ENV{LAPACKE_LIBDIR}")
    if(ENV_LAPACKE_LIBDIR)
      list(APPEND _lib_env "${ENV_LAPACKE_LIBDIR}")
    elseif(ENV_LAPACKE_DIR)
      list(APPEND _lib_env "${ENV_LAPACKE_DIR}")
      list(APPEND _lib_env "${ENV_LAPACKE_DIR}/lib")
    else()
      if(WIN32)
        string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
      else()
        if(APPLE)
          string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
        else()
          string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
        endif()
        list(APPEND _lib_env "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
        list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
      endif()
    endif()
    list(REMOVE_DUPLICATES _lib_env)

    # Try to find the lapacke lib in the given paths
    # ----------------------------------------------

    # name of the lapacke library
    set(LAPACKE_lapacke_NAMES "lapacke")
    if(LAPACKE_STATIC)
        if(WIN32)
            set(LAPACKE_lapacke_NAMES "liblapacke.lib")
        endif()

        if(UNIX)
            set(LAPACKE_lapacke_NAMES "liblapacke.a")
        endif()
    endif()

    # call cmake macro to find the lib path
    if(LAPACKE_LIBDIR)
      set(LAPACKE_lapacke_LIBRARY "LAPACKE_lapacke_LIBRARY-NOTFOUND")
      find_library(LAPACKE_lapacke_LIBRARY
        NAMES ${LAPACKE_lapacke_NAMES}
        HINTS ${LAPACKE_LIBDIR})
    else()
      if(LAPACKE_DIR)
        set(LAPACKE_lapacke_LIBRARY "LAPACKE_lapacke_LIBRARY-NOTFOUND")
        find_library(LAPACKE_lapacke_LIBRARY
          NAMES ${LAPACKE_lapacke_NAMES}
          HINTS ${LAPACKE_DIR}
          PATH_SUFFIXES lib lib32 lib64)
      else()
        set(LAPACKE_lapacke_LIBRARY "LAPACKE_lapacke_LIBRARY-NOTFOUND")
        find_library(LAPACKE_lapacke_LIBRARY
          NAMES ${LAPACKE_lapacke_NAMES}
          HINTS ${_lib_env})
      endif()
    endif()
    mark_as_advanced(LAPACKE_lapacke_LIBRARY)

    # If found, add path to cmake variable
    # ------------------------------------
    if (LAPACKE_lapacke_LIBRARY)
      get_filename_component(lapacke_lib_path "${LAPACKE_lapacke_LIBRARY}" PATH)
      # set cmake variables
      set(LAPACKE_LIBRARIES    "${LAPACKE_lapacke_LIBRARY}")
      set(LAPACKE_LIBRARY_DIRS "${lapacke_lib_path}")
    else ()
      set(LAPACKE_LIBRARIES    "LAPACKE_LIBRARIES-NOTFOUND")
      set(LAPACKE_LIBRARY_DIRS "LAPACKE_LIBRARY_DIRS-NOTFOUND")
      if (NOT LAPACKE_FIND_QUIETLY)
        message(STATUS "Looking for lapacke -- lib lapacke not found")
      endif()
    endif ()

    # check a function to validate the find
    if(LAPACKE_LIBRARIES)

      set(REQUIRED_LDFLAGS)
      set(REQUIRED_INCDIRS)
      set(REQUIRED_LIBDIRS)
      set(REQUIRED_LIBS)

      # LAPACKE
      if (LAPACKE_INCLUDE_DIRS)
        set(REQUIRED_INCDIRS "${LAPACKE_INCLUDE_DIRS}")
      endif()
      if (LAPACKE_LIBRARY_DIRS)
        set(REQUIRED_LIBDIRS "${LAPACKE_LIBRARY_DIRS}")
      endif()
      set(REQUIRED_LIBS "${LAPACKE_LIBRARIES}")
      # LAPACK
      if (LAPACK_INCLUDE_DIRS)
        list(APPEND REQUIRED_INCDIRS "${LAPACK_INCLUDE_DIRS}")
      endif()
      if (LAPACK_LIBRARY_DIRS)
        list(APPEND REQUIRED_LIBDIRS "${LAPACK_LIBRARY_DIRS}")
      endif()
      list(APPEND REQUIRED_LIBS "${LAPACK_LIBRARIES}")
      if (LAPACK_LINKER_FLAGS)
        list(APPEND REQUIRED_LDFLAGS "${LAPACK_LINKER_FLAGS}")
      endif()
      # Fortran
      if (CMAKE_C_COMPILER_ID MATCHES "GNU")
        find_library(
          FORTRAN_gfortran_LIBRARY
          NAMES gfortran
          HINTS ${_lib_env}
          )
        mark_as_advanced(FORTRAN_gfortran_LIBRARY)
        if (FORTRAN_gfortran_LIBRARY)
          list(APPEND REQUIRED_LIBS "${FORTRAN_gfortran_LIBRARY}")
        endif()
      elseif (CMAKE_C_COMPILER_ID MATCHES "Intel")
        find_library(
          FORTRAN_ifcore_LIBRARY
          NAMES ifcore
          HINTS ${_lib_env}
          )
        mark_as_advanced(FORTRAN_ifcore_LIBRARY)
        if (FORTRAN_ifcore_LIBRARY)
          list(APPEND REQUIRED_LIBS "${FORTRAN_ifcore_LIBRARY}")
        endif()
      endif()
      # m
      find_library(M_LIBRARY NAMES m HINTS ${_lib_env})
      mark_as_advanced(M_LIBRARY)
      if(M_LIBRARY)
        list(APPEND REQUIRED_LIBS "-lm")
      endif()
      # set required libraries for link
      set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
      set(CMAKE_REQUIRED_LIBRARIES)
      list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LDFLAGS}")
      foreach(lib_dir ${REQUIRED_LIBDIRS})
        list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
      endforeach()
      list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
      string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

      # test link
      unset(LAPACKE_WORKS CACHE)
      include(CheckFunctionExists)
      check_function_exists(LAPACKE_dgeqrf LAPACKE_WORKS)
      mark_as_advanced(LAPACKE_WORKS)

      if(LAPACKE_WORKS)
        # save link with dependencies
        set(LAPACKE_LIBRARIES_DEP "${REQUIRED_LIBS}")
        set(LAPACKE_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
        set(LAPACKE_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
        set(LAPACKE_LINKER_FLAGS "${REQUIRED_LDFLAGS}")
        list(REMOVE_DUPLICATES LAPACKE_LIBRARY_DIRS_DEP)
        list(REMOVE_DUPLICATES LAPACKE_INCLUDE_DIRS_DEP)
        list(REMOVE_DUPLICATES LAPACKE_LINKER_FLAGS)
      else()
        if(NOT LAPACKE_FIND_QUIETLY)
          message(STATUS "Looking for lapacke: test of LAPACKE_dgeqrf with lapacke and lapack libraries fails")
          message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
          message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
          message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
        endif()
      endif()
      set(CMAKE_REQUIRED_INCLUDES)
      set(CMAKE_REQUIRED_FLAGS)
      set(CMAKE_REQUIRED_LIBRARIES)
    endif(LAPACKE_LIBRARIES)


else(ACCELERATE)

  if (NOT LAPACKE_FIND_QUIETLY)
    message(STATUS "Eigen Overload requires Apple Accelerate and it has not been found.")
  endif()

endif(ACCELERATE)

if (LAPACKE_LIBRARIES)
  list(GET LAPACKE_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (${first_lib_path} MATCHES "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)")
    string(REGEX REPLACE "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)" "" not_cached_dir "${first_lib_path}")
    set(LAPACKE_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of LAPACKE library" FORCE)
  else()
    set(LAPACKE_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of LAPACKE library" FORCE)
  endif()
endif()
mark_as_advanced(LAPACKE_DIR)
mark_as_advanced(LAPACKE_DIR_FOUND)

list(PREPEND LAPACKE_LIBRARIES "${ACCELERATE}")

# check that APPLE_ACCELERATE has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(APPLE_ACCELERATE DEFAULT_MSG
  LAPACKE_LIBRARIES
  LAPACKE_WORKS)


  # ---------------------------------------------------------------------------
# Imported targets
# ---------------------------------------------------------------------------
if(APPLE_ACCELERATE_FOUND)

  # 1) Accelerate framework target
  if(ACCELERATE AND NOT TARGET APPLE_ACCELERATE::Accelerate)
    add_library(APPLE_ACCELERATE::Accelerate INTERFACE IMPORTED)
    set_target_properties(APPLE_ACCELERATE::Accelerate PROPERTIES
      INTERFACE_LINK_LIBRARIES "${ACCELERATE}"
    )
  endif()

  # 2) LAPACKE (+deps) target, for Eigen overload usage
  #
  # Prefer the "validated" *_DEP variables when available (they include deps).
  if(LAPACKE_WORKS AND LAPACKE_LIBRARIES AND NOT TARGET APPLE_ACCELERATE::LAPACKE)

    add_library(APPLE_ACCELERATE::LAPACKE INTERFACE IMPORTED)

    # Includes
    if(DEFINED LAPACKE_INCLUDE_DIRS_DEP AND LAPACKE_INCLUDE_DIRS_DEP)
      set_property(TARGET APPLE_ACCELERATE::LAPACKE PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES "${LAPACKE_INCLUDE_DIRS_DEP}"
      )
    elseif(DEFINED LAPACKE_INCLUDE_DIRS AND LAPACKE_INCLUDE_DIRS)
      set_property(TARGET APPLE_ACCELERATE::LAPACKE PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES "${LAPACKE_INCLUDE_DIRS}"
      )
    endif()

    # Link libs (use the “DEP” set if present, otherwise fall back)
    if(DEFINED LAPACKE_LIBRARIES_DEP AND LAPACKE_LIBRARIES_DEP)
      set_property(TARGET APPLE_ACCELERATE::LAPACKE PROPERTY
        INTERFACE_LINK_LIBRARIES "${LAPACKE_LIBRARIES_DEP}"
      )
    else()
      set_property(TARGET APPLE_ACCELERATE::LAPACKE PROPERTY
        INTERFACE_LINK_LIBRARIES "${LAPACKE_LIBRARIES}"
      )
    endif()

    # Link options / flags, if any
    if(DEFINED LAPACKE_LINKER_FLAGS AND LAPACKE_LINKER_FLAGS)
      # CMake >= 3.13 supports INTERFACE_LINK_OPTIONS
      set_property(TARGET APPLE_ACCELERATE::LAPACKE APPEND PROPERTY
        INTERFACE_LINK_OPTIONS "${LAPACKE_LINKER_FLAGS}"
      )
    endif()

    # Ensure Accelerate is in the chain (your module prepends it to LAPACKE_LIBRARIES,
    # but in case you end up using *_DEP which might not include it)
    if(TARGET APPLE_ACCELERATE::Accelerate)
      set_property(TARGET APPLE_ACCELERATE::LAPACKE APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES APPLE_ACCELERATE::Accelerate
      )
    endif()

  endif()

endif()
