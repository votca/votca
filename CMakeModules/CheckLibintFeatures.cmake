# Reports whether the libint2 found for XTP has the derivative integrals that
# the DFT engine's compute_forces option needs: 1st derivatives of the
# one-body, 2-centre and 3-centre integrals (libint configured with
# --enable-1body=1 --enable-eri2=1 --enable-eri3=1, or higher).
#
# Reads libint2/config.h from the include directories of Libint2::int2; no
# try_compile. libint <= 2.9 names the macros INCLUDE_ONEBODY/_ERI2/_ERI3,
# later versions LIBINT_INCLUDE_ONEBODY/_ERI2/_ERI3; both are recognised.
#
# Sets VOTCA_LIBINT_FORCES (ON/OFF), VOTCA_LIBINT_CONFIG_H (also cached, for
# reference) and VOTCA_LIBINT_EXPECT_DEFS (compile definitions for the
# header check in libint2_derivative_calls.cc).
# With -DREQUIRE_LIBINT_FORCES=ON, missing derivatives are a configure error.

option(REQUIRE_LIBINT_FORCES "Fail at configure time if libint2 lacks the derivative integrals for analytic DFT forces" OFF)

function(votca_check_libint_features)
  set(_cfg "")
  get_target_property(_incs Libint2::int2 INTERFACE_INCLUDE_DIRECTORIES)
  if(NOT _incs)
    set(_incs "")
  endif()
  # libint's own config packages put the C++ headers on a separate target
  foreach(_t Libint2::cxx Libint2::int2-cxx)
    if(TARGET ${_t})
      get_target_property(_more ${_t} INTERFACE_INCLUDE_DIRECTORIES)
      if(_more)
        list(APPEND _incs ${_more})
      endif()
    endif()
  endforeach()
  foreach(_d IN LISTS _incs)
    string(REGEX REPLACE "^\\$<BUILD_INTERFACE:(.*)>$" "\\1" _d "${_d}")
    string(REGEX REPLACE "^\\$<INSTALL_INTERFACE:.*>$" "" _d "${_d}")
    if(_d AND EXISTS "${_d}/libint2/config.h")
      set(_cfg "${_d}/libint2/config.h")
      break()
    endif()
  endforeach()

  # BUILD_OWN_LIBINT before its first build: read config.h from the extracted export
  if(NOT _cfg AND VOTCA_LIBINT_PENDING_CONFIG_H AND EXISTS "${VOTCA_LIBINT_PENDING_CONFIG_H}")
    set(_cfg "${VOTCA_LIBINT_PENDING_CONFIG_H}")
  endif()

  set(VOTCA_LIBINT_CONFIG_H "${_cfg}" PARENT_SCOPE)
  set(VOTCA_LIBINT_CONFIG_H "${_cfg}" CACHE INTERNAL "libint2/config.h of the libint2 used by XTP")
  if(NOT _cfg)
    set(VOTCA_LIBINT_FORCES OFF PARENT_SCOPE)
    message(STATUS "libint2 derivative integrals: unknown (libint2/config.h not found in ${_incs}); "
                   "analytic DFT forces are checked when compiling libint2_derivative_calls.cc")
    if(REQUIRE_LIBINT_FORCES)
      message(FATAL_ERROR "REQUIRE_LIBINT_FORCES=ON, but libint2/config.h was not found")
    endif()
    return()
  endif()

  file(STRINGS "${_cfg}" _lines REGEX "^#define[ \t]+(LIBINT_)?INCLUDE_(ONEBODY|ERI|ERI2|ERI3)[ \t]")
  set(_missing "")
  set(_summary "")
  set(_expect_orders "")
  foreach(_k ONEBODY ERI2 ERI3)
    set(_order -1)
    foreach(_l IN LISTS _lines)
      if(_l MATCHES "^#define[ \t]+(LIBINT_)?INCLUDE_${_k}[ \t]+(-?[0-9]+)")
        set(_order ${CMAKE_MATCH_2})
      endif()
    endforeach()
    list(APPEND _summary "${_k}=${_order}")
    list(APPEND _expect_orders "VOTCA_LIBINT_EXPECT_${_k}=${_order}")
    if(_order LESS 1)
      list(APPEND _missing ${_k})
    endif()
  endforeach()
  string(REPLACE ";" ", " _summary "${_summary}")
  # 4-centre ERIs: not needed for forces (XTP uses RI), but their derivative
  # order decides which libint functions the headers reference, so it is part
  # of the header check below
  set(_eri4 -1)
  foreach(_l IN LISTS _lines)
    if(_l MATCHES "^#define[ \t]+(LIBINT_)?INCLUDE_ERI[ \t]+(-?[0-9]+)")
      set(_eri4 ${CMAKE_MATCH_2})
    endif()
  endforeach()
  list(APPEND _expect_orders "VOTCA_LIBINT_EXPECT_ERI=${_eri4}")

  # version and derivative orders, for the compile-time check in
  # libint2_derivative_calls.cc that the compiler sees the same headers
  file(STRINGS "${_cfg}" _vlines REGEX "^#define[ \t]+LIBINT_(MAJOR|MINOR|MICRO)_VERSION[ \t]")
  set(_defs "")
  foreach(_p MAJOR MINOR MICRO)
    foreach(_l IN LISTS _vlines)
      if(_l MATCHES "^#define[ \t]+LIBINT_${_p}_VERSION[ \t]+([0-9]+)")
        list(APPEND _defs "VOTCA_LIBINT_EXPECT_${_p}=${CMAKE_MATCH_1}")
      endif()
    endforeach()
  endforeach()
  list(LENGTH _defs _nv)
  if(_nv EQUAL 3)
    list(APPEND _defs ${_expect_orders})
    set(VOTCA_LIBINT_EXPECT_DEFS "${_defs}" PARENT_SCOPE)
  endif()

  if(APPLE)
    set(_hint_cache "~/Library/Caches/votca")
  else()
    set(_hint_cache "~/.cache/votca")
  endif()
  if(_missing)
    set(VOTCA_LIBINT_FORCES OFF PARENT_SCOPE)
    if(BUILD_OWN_LIBINT)
      set(_hint "the export in VOTCA_LIBINT_EXPORT_URL lacks them; leave it empty for the VOTCA default export")
    else()
      set(_hint "-DBUILD_OWN_LIBINT=ON builds a libint2 with them")
    endif()
    message(STATUS "libint2 derivative integrals: ${_summary} (${_cfg}) -> analytic DFT forces OFF (${_hint})")
    if(REQUIRE_LIBINT_FORCES)
      if(BUILD_OWN_LIBINT)
        message(FATAL_ERROR
          "REQUIRE_LIBINT_FORCES=ON, but the libint2 export ${VOTCA_LIBINT_EXPORT_URL} has no "
          "1st derivatives for: ${_missing}. Leave VOTCA_LIBINT_EXPORT_URL empty to use the VOTCA "
          "default export, or generate one with CMakeModules/make_libint_export.sh.")
      else()
        message(FATAL_ERROR
          "REQUIRE_LIBINT_FORCES=ON, but the libint2 found (${_cfg}) has no 1st derivatives for: "
          "${_missing}. Configure with -DBUILD_OWN_LIBINT=ON instead: it downloads the VOTCA libint "
          "export (with derivatives) and builds it once per machine into ${_hint_cache}. "
          "Alternatively, point CMAKE_PREFIX_PATH to a libint2 configured with "
          "LIBINT2_ENABLE_ONEBODY=1, LIBINT2_ENABLE_ERI2=1 and LIBINT2_ENABLE_ERI3=1.")
      endif()
    endif()
  else()
    set(VOTCA_LIBINT_FORCES ON PARENT_SCOPE)
    message(STATUS "libint2 derivative integrals: ${_summary} (${_cfg}) -> analytic DFT forces ON")
  endif()
endfunction()
