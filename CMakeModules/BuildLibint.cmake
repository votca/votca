# Builds libint2 from a pre-generated library export, once per machine.
#
# A libint "export" is the generated integral library (libint-<ver>.tgz,
# produced by the libint generator with `make export`). It builds with CMake
# in minutes and needs no autotools, GMP or generator run. XTP needs an export
# generated with derivative integrals (--enable-1body=1 --enable-eri2=1
# --enable-eri3=1); the official release exports do not have them. See
# CMakeModules/make_libint_export.sh for how the VOTCA export is generated.
#
# Configure only downloads (SHA256-checked) and extracts the export into a
# per-machine cache outside the build tree:
#   ${VOTCA_LIBINT_CACHE_DIR}/libint-<export sha256 prefix>-<compiler>-<arch>
# and defines Libint2::int2 at its install location there. The library itself
# is compiled by the first `cmake --build` (target votca_libint, which runs
# before the rest of VOTCA, so the outputs do not interleave). Once it is
# installed, every later configure, in any build directory and on any VOTCA
# branch, uses the installed library and adds no build step. The cache key contains only the
# export's checksum, the compiler and the architecture, so changing VOTCA
# sources or switching branches never rebuilds libint. A file lock serialises
# concurrent configures and builds; a stamp file marks a complete install, so
# an interrupted build is redone by the next build.
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

set(_srcroot "${VOTCA_LIBINT_CACHE_DIR}/${_key}.src")
set(_lock "${VOTCA_LIBINT_CACHE_DIR}/${_key}.lock")

if(EXISTS "${_stamp}")
  # ---- installed: use libint's own CMake package ---------------------------
  message(STATUS "libint2: using the cached build in ${VOTCA_LIBINT_PREFIX}")
  find_package(Libint2 CONFIG REQUIRED PATHS "${VOTCA_LIBINT_PREFIX}" NO_DEFAULT_PATH)
  if(NOT TARGET Libint2::int2)
    message(FATAL_ERROR "libint2 in ${VOTCA_LIBINT_PREFIX} provides no Libint2::int2 target")
  endif()
  return()
endif()

# ---- not installed yet: download + extract now (seconds), compile at build time
file(MAKE_DIRECTORY "${VOTCA_LIBINT_CACHE_DIR}")
file(LOCK "${_lock}" GUARD FILE TIMEOUT 7200 RESULT_VARIABLE _lockres)
if(NOT _lockres EQUAL 0)
  message(FATAL_ERROR "Could not lock ${_lock}: ${_lockres}")
endif()
if(NOT EXISTS "${_srcroot}/.extracted")
  file(REMOVE_RECURSE "${_srcroot}")
  file(MAKE_DIRECTORY "${_srcroot}")
  if(NOT _tarball)
    set(_tarball "${_srcroot}/libint-export.tgz")
    message(STATUS "libint2: downloading ${VOTCA_LIBINT_EXPORT_URL} (once per machine)")
    file(DOWNLOAD "${VOTCA_LIBINT_EXPORT_URL}" "${_tarball}"
      EXPECTED_HASH SHA256=${_sha} TLS_VERIFY ON SHOW_PROGRESS STATUS _dl)
    list(GET _dl 0 _dlcode)
    if(NOT _dlcode EQUAL 0)
      message(FATAL_ERROR "libint2: download failed: ${_dl}")
    endif()
  endif()
  file(ARCHIVE_EXTRACT INPUT "${_tarball}" DESTINATION "${_srcroot}")
  if(_tarball STREQUAL "${_srcroot}/libint-export.tgz")
    file(REMOVE "${_tarball}")
  endif()
  # exports made on macOS can carry AppleDouble metadata files (._*); drop them
  file(GLOB_RECURSE _appledouble "${_srcroot}/._*" "${_srcroot}/*/._*")
  if(_appledouble)
    file(REMOVE ${_appledouble})
  endif()
  file(GLOB _dirs LIST_DIRECTORIES true "${_srcroot}/*")
  list(LENGTH _dirs _n)
  if(NOT _n EQUAL 1 OR NOT EXISTS "${_dirs}/CMakeLists.txt")
    message(FATAL_ERROR "libint2: ${VOTCA_LIBINT_EXPORT_URL} is not a libint export "
                        "(expected one top-level directory with a CMakeLists.txt)")
  endif()
  file(WRITE "${_srcroot}/.extracted" "${_dirs}\n")
endif()
file(LOCK "${_lock}" RELEASE)
file(STRINGS "${_srcroot}/.extracted" _srcdir LIMIT_COUNT 1)

# the export ships its generated config.h; CheckLibintFeatures reads it from
# here because the install location is still empty at configure time
set(VOTCA_LIBINT_PENDING_CONFIG_H "${_srcdir}/include/libint2/config.h")

if(VOTCA_LIBINT_BUILD_JOBS)
  set(_jobs ${VOTCA_LIBINT_BUILD_JOBS})
else()
  cmake_host_system_information(RESULT _jobs QUERY NUMBER_OF_LOGICAL_CORES)
endif()

# settings for the libint build, passed through an initial-cache file (keeps lists intact)
set(_init "")
foreach(_v CMAKE_CXX_COMPILER CMAKE_C_COMPILER CMAKE_PREFIX_PATH Eigen3_DIR Boost_DIR
           CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT)
  if(${_v})
    string(APPEND _init "set(${_v} \"${${_v}}\" CACHE STRING \"\")\n")
  endif()
endforeach()
file(WRITE "${CMAKE_BINARY_DIR}/libint/init.cmake" "${_init}")

message(STATUS "libint2: will be compiled by the first build (target votca_libint) and installed into "
               "${VOTCA_LIBINT_PREFIX}; later configures reuse it")
add_custom_command(
  OUTPUT "${_stamp}"
  BYPRODUCTS "${VOTCA_LIBINT_PREFIX}/lib/libint2.a"
  COMMAND ${CMAKE_COMMAND}
    "-DSRC=${_srcdir}" "-DPREFIX=${VOTCA_LIBINT_PREFIX}" "-DSTAMP=${_stamp}" "-DLOCK=${_lock}"
    "-DWORK=${VOTCA_LIBINT_CACHE_DIR}/${_key}.build" "-DINIT=${CMAKE_BINARY_DIR}/libint/init.cmake"
    "-DGENERATOR=${CMAKE_GENERATOR}" "-DJOBS=${_jobs}" "-DSHA=${_sha}"
    -P "${CMAKE_CURRENT_LIST_DIR}/BuildLibintStep.cmake"
  DEPENDS "${CMAKE_CURRENT_LIST_DIR}/BuildLibintStep.cmake"
  USES_TERMINAL
  COMMENT "libint2: one-time build for this machine (${_jobs} jobs); VOTCA builds afterwards")
add_custom_target(votca_libint DEPENDS "${_stamp}")

# Libint2::int2 at its future install location, with the properties libint's
# own package gives Libint2::int-static
file(MAKE_DIRECTORY "${VOTCA_LIBINT_PREFIX}/include")
add_library(Libint2::int2 STATIC IMPORTED)
set_target_properties(Libint2::int2 PROPERTIES
  IMPORTED_LOCATION "${VOTCA_LIBINT_PREFIX}/lib/libint2.a"
  IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
  INTERFACE_INCLUDE_DIRECTORIES "${VOTCA_LIBINT_PREFIX}/include"
  INTERFACE_COMPILE_FEATURES "cxx_std_11")
add_dependencies(Libint2::int2 votca_libint)
