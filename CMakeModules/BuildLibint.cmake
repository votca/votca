# Builds libint2 from a pre-generated library export, once per machine.
#
# A libint "export" is the generated integral library (libint-<ver>.tgz,
# produced by the libint generator with `make export`). It builds with CMake
# in minutes and needs no autotools, GMP or generator run. XTP needs an export
# generated with derivative integrals (--enable-1body=1 --enable-eri2=1
# --enable-eri3=1); the official release exports do not have them. See
# CMakeModules/make_libint_export.sh for how the VOTCA export is generated.
#
# The library is built at configure time and installed into a per-machine
# cache outside the build tree:
#   ${VOTCA_LIBINT_CACHE_DIR}/libint-<export sha256 prefix>-<compiler>-<arch>
# Every later configure, in any build directory and on any VOTCA branch,
# finds the installed library there and builds nothing. The cache key
# contains only the export's checksum, the compiler and the architecture, so
# changing VOTCA sources or switching branches never rebuilds libint.
# A file lock serialises concurrent configures; a stamp file marks a complete
# install, so an interrupted build is redone on the next configure.
#
# Variables:
#   VOTCA_LIBINT_EXPORT_URL     URL or local path of the export tarball (default: VOTCA release asset)
#   VOTCA_LIBINT_EXPORT_SHA256  its SHA256 (required for URLs, it is the cache key)
#   VOTCA_LIBINT_CACHE_DIR      default: $VOTCA_LIBINT_CACHE_DIR, else
#                               $XDG_CACHE_HOME/votca, else ~/.cache/votca
#                               (~/Library/Caches/votca on macOS)
#   VOTCA_LIBINT_BUILD_JOBS     parallel jobs for the one-time build (default: all cores)
#
# Remove the cache directory to force a rebuild.

# Default: libint 2.13.1 with 1st derivatives of the 1-body, 2-centre and
# 3-centre integrals, generated with make_libint_export.sh (128 MB).
set(VOTCA_LIBINT_DEFAULT_URL
  "https://github.com/votca/libint-export/releases/download/v2.13.1-deriv1/libint-2.13.1-votca-deriv1.tgz")
set(VOTCA_LIBINT_DEFAULT_SHA256 "468683b3575f050698391a5d4d5607f3e728353c4153b710939289e215a89bd2")

set(VOTCA_LIBINT_EXPORT_URL "" CACHE STRING
  "URL or path of a pre-generated libint2 export tarball (empty: the VOTCA default export)")
set(VOTCA_LIBINT_EXPORT_SHA256 "" CACHE STRING "SHA256 of the libint2 export tarball")

if(NOT VOTCA_LIBINT_EXPORT_URL)
  # (-DVOTCA_LIBINT_EXPORT_URL=<local copy of the tarball> avoids the download)
  set(VOTCA_LIBINT_EXPORT_URL "${VOTCA_LIBINT_DEFAULT_URL}")
  set(VOTCA_LIBINT_EXPORT_SHA256 "${VOTCA_LIBINT_DEFAULT_SHA256}")
endif()

# ---- cache location -------------------------------------------------------
if(NOT VOTCA_LIBINT_CACHE_DIR)
  if(DEFINED ENV{VOTCA_LIBINT_CACHE_DIR})
    set(_cache "$ENV{VOTCA_LIBINT_CACHE_DIR}")
  elseif(DEFINED ENV{XDG_CACHE_HOME})
    set(_cache "$ENV{XDG_CACHE_HOME}/votca")
  elseif(APPLE)
    set(_cache "$ENV{HOME}/Library/Caches/votca")
  else()
    set(_cache "$ENV{HOME}/.cache/votca")
  endif()
  set(VOTCA_LIBINT_CACHE_DIR "${_cache}" CACHE PATH "Per-machine cache for the libint2 build")
endif()

# ---- source: local file or download ---------------------------------------
if(EXISTS "${VOTCA_LIBINT_EXPORT_URL}")
  get_filename_component(_tarball "${VOTCA_LIBINT_EXPORT_URL}" ABSOLUTE)
  file(SHA256 "${_tarball}" _sha)
  if(VOTCA_LIBINT_EXPORT_SHA256 AND NOT _sha STREQUAL VOTCA_LIBINT_EXPORT_SHA256)
    message(FATAL_ERROR "SHA256 of ${_tarball} is ${_sha}, expected ${VOTCA_LIBINT_EXPORT_SHA256}")
  endif()
else()
  if(NOT VOTCA_LIBINT_EXPORT_SHA256)
    message(FATAL_ERROR "VOTCA_LIBINT_EXPORT_SHA256 is required when VOTCA_LIBINT_EXPORT_URL is a URL")
  endif()
  set(_sha "${VOTCA_LIBINT_EXPORT_SHA256}")
  set(_tarball "")
endif()

# ---- cache key: export checksum + compiler + architecture -----------------
string(SUBSTRING "${_sha}" 0 12 _sha12)
string(REGEX MATCH "^[0-9]+" _cxx_major "${CMAKE_CXX_COMPILER_VERSION}")
set(_arch "${CMAKE_SYSTEM_PROCESSOR}")
if(APPLE AND CMAKE_OSX_ARCHITECTURES)
  string(REPLACE ";" "+" _arch "${CMAKE_OSX_ARCHITECTURES}")
endif()
string(TOLOWER "libint-${_sha12}-${CMAKE_CXX_COMPILER_ID}${_cxx_major}-${CMAKE_SYSTEM_NAME}-${_arch}" _key)
set(VOTCA_LIBINT_PREFIX "${VOTCA_LIBINT_CACHE_DIR}/${_key}")
set(_stamp "${VOTCA_LIBINT_PREFIX}/.votca-libint-complete")

# ---- build once -------------------------------------------------------------
if(NOT EXISTS "${_stamp}")
  file(MAKE_DIRECTORY "${VOTCA_LIBINT_CACHE_DIR}")
  file(LOCK "${VOTCA_LIBINT_CACHE_DIR}/${_key}.lock" GUARD FILE TIMEOUT 7200 RESULT_VARIABLE _lockres)
  if(NOT _lockres EQUAL 0)
    message(FATAL_ERROR "Could not lock ${VOTCA_LIBINT_CACHE_DIR}/${_key}.lock: ${_lockres}")
  endif()
  # another configure may have finished the build while we waited for the lock
  if(NOT EXISTS "${_stamp}")
    set(_work "${VOTCA_LIBINT_CACHE_DIR}/${_key}.work")
    file(REMOVE_RECURSE "${_work}" "${VOTCA_LIBINT_PREFIX}")
    file(MAKE_DIRECTORY "${_work}/src")
    if(NOT _tarball)
      set(_tarball "${_work}/libint-export.tgz")
      message(STATUS "libint2: downloading ${VOTCA_LIBINT_EXPORT_URL}")
      file(DOWNLOAD "${VOTCA_LIBINT_EXPORT_URL}" "${_tarball}"
        EXPECTED_HASH SHA256=${_sha} TLS_VERIFY ON STATUS _dl)
      list(GET _dl 0 _dlcode)
      if(NOT _dlcode EQUAL 0)
        message(FATAL_ERROR "libint2: download failed: ${_dl}")
      endif()
    endif()
    file(ARCHIVE_EXTRACT INPUT "${_tarball}" DESTINATION "${_work}/src")
    # exports made on macOS can carry AppleDouble metadata files (._*); drop them
    file(GLOB_RECURSE _appledouble "${_work}/src/._*" "${_work}/src/*/._*")
    if(_appledouble)
      file(REMOVE ${_appledouble})
    endif()
    file(GLOB _srcdir LIST_DIRECTORIES true "${_work}/src/*")
    list(LENGTH _srcdir _n)
    if(NOT _n EQUAL 1 OR NOT EXISTS "${_srcdir}/CMakeLists.txt")
      message(FATAL_ERROR "libint2: ${_tarball} is not a libint export (expected one top-level directory with a CMakeLists.txt)")
    endif()

    if(VOTCA_LIBINT_BUILD_JOBS)
      set(_jobs ${VOTCA_LIBINT_BUILD_JOBS})
    else()
      cmake_host_system_information(RESULT _jobs QUERY NUMBER_OF_LOGICAL_CORES)
    endif()
    message(STATUS "libint2: one-time build into ${VOTCA_LIBINT_PREFIX} (${_jobs} jobs; "
                   "log in ${_work}/build.log). Later configures reuse it.")
    set(_args
      -S "${_srcdir}" -B "${_work}/build"
      -G "${CMAKE_GENERATOR}"
      -DCMAKE_BUILD_TYPE=Release
      -DCMAKE_INSTALL_PREFIX=${VOTCA_LIBINT_PREFIX}
      -DCMAKE_INSTALL_LIBDIR=lib
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DBUILD_SHARED_LIBS=OFF
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON
      -DLIBINT2_REQUIRE_CXX_API=ON
      -DLIBINT2_REQUIRE_CXX_API_COMPILED=OFF
      -DLIBINT2_ENABLE_FORTRAN=OFF
      -DLIBINT2_ENABLE_PYTHON=OFF
      -DBUILD_TESTING=OFF)
    # list-valued settings go through an initial-cache file, which keeps the ';'
    set(_init "")
    foreach(_v CMAKE_C_COMPILER CMAKE_PREFIX_PATH Eigen3_DIR Boost_DIR CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET)
      if(${_v})
        string(APPEND _init "set(${_v} \"${${_v}}\" CACHE STRING \"\")\n")
      endif()
    endforeach()
    file(WRITE "${_work}/init.cmake" "${_init}")
    list(APPEND _args -C "${_work}/init.cmake")
    function(_votca_libint_step)
      execute_process(COMMAND ${ARGN} RESULT_VARIABLE _rc
        OUTPUT_FILE "${_work}/step.log" ERROR_FILE "${_work}/step.log")
      file(READ "${_work}/step.log" _log)
      file(APPEND "${_work}/build.log" "${_log}")
      if(NOT _rc EQUAL 0)
        message(FATAL_ERROR "libint2: build step failed (${_rc}); see ${_work}/build.log")
      endif()
    endfunction()
    _votca_libint_step(${CMAKE_COMMAND} ${_args})
    _votca_libint_step(${CMAKE_COMMAND} --build "${_work}/build" --parallel ${_jobs} --config Release)
    _votca_libint_step(${CMAKE_COMMAND} --install "${_work}/build" --config Release)
    file(WRITE "${_stamp}" "export=${VOTCA_LIBINT_EXPORT_URL}\nsha256=${_sha}\ncompiler=${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_VERSION}\n")
    file(REMOVE_RECURSE "${_work}")
  else()
    message(STATUS "libint2: another configure built ${VOTCA_LIBINT_PREFIX}; using it")
  endif()
  file(LOCK "${VOTCA_LIBINT_CACHE_DIR}/${_key}.lock" RELEASE)
else()
  message(STATUS "libint2: using the cached build in ${VOTCA_LIBINT_PREFIX}")
endif()

find_package(Libint2 CONFIG REQUIRED PATHS "${VOTCA_LIBINT_PREFIX}" NO_DEFAULT_PATH)
if(NOT TARGET Libint2::int2)
  message(FATAL_ERROR "libint2 in ${VOTCA_LIBINT_PREFIX} provides no Libint2::int2 target")
endif()
