########################################################################
# User input options                                                   #
########################################################################
find_package(Libxc QUIET CONFIG)
set_package_properties(Libxc PROPERTIES TYPE REQUIRED PURPOSE "Provides exchange-correclation functionals")
if(NOT Libxc_FOUND)
  find_package(Libxc MODULE REQUIRED)
endif()
if(NOT TARGET Libint2::int2)
  find_package(libint2 QUIET CONFIG)
  set_package_properties(libint2 PROPERTIES TYPE REQUIRED PURPOSE "Calculates Gaussian integrals")
  if(NOT libint2_FOUND)
    find_package(libint2 2.6 MODULE REQUIRED)
  endif()
endif()

# Report, at configure time, whether this libint2 was built with
# derivative-integral support for each of the three independent
# operator categories libint2_derivative_calls.cc actually uses --
# LIBINT_INCLUDE_ONEBODY (overlap, kinetic, nuclear attraction),
# LIBINT_INCLUDE_ERI2 (the two-center Coulomb metric),
# LIBINT_INCLUDE_ERI3 (three-center RI) -- each set independently at
# libint2's OWN configure time (LIBINT2_ENABLE_ONEBODY/_ERI2/_ERI3).
# Deliberately check_cxx_symbol_exists, NOT try_run: this only needs to
# COMPILE a one-line check (does #include <libint2.hpp> define this
# macro), never actually running anything -- avoiding the exact
# category of fragility try_run had (a full, standalone test program
# with its own bugs to get right: an earlier version of this check hit
# both a libint2 API mismatch and a missing Eigen include path before
# being replaced with the current, simpler compile-time #ifdef guards
# in libint2_derivative_calls.cc itself; this CMake-level check exists
# only to surface the same information as a clear, early message,
# purely for visibility -- the actual C++ guards do not depend on it).
include(CheckCXXSymbolExists)
# CMAKE_REQUIRED_LIBRARIES accepts imported targets directly (modern
# CMake) and pulls in their FULL interface -- include directories,
# compile definitions, and transitive dependencies -- automatically.
# Deliberately NOT get_target_property(... INTERFACE_INCLUDE_DIRECTORIES)
# on each target individually: that approach was tried first and
# produced a false "not found" for all three symbols even on a machine
# with confirmed-working derivative support, most likely because
# libint2's own headers transitively need Eigen/Core (confirmed
# separately: this exact "'Eigen/Core' file not found" error showed up
# during the earlier, abandoned try_run attempt) and manually extracted
# include-directory properties do not reliably capture everything a
# full target-based CMAKE_REQUIRED_LIBRARIES resolution does.
set(CMAKE_REQUIRED_LIBRARIES Libint2::int2 Eigen3::Eigen)
set(CMAKE_REQUIRED_FLAGS "-std=c++17")
check_cxx_symbol_exists(LIBINT_INCLUDE_ONEBODY "libint2.hpp" XTP_LIBINT_HAS_ONEBODY_DERIVS)
check_cxx_symbol_exists(LIBINT_INCLUDE_ERI2 "libint2.hpp" XTP_LIBINT_HAS_ERI2_DERIVS)
check_cxx_symbol_exists(LIBINT_INCLUDE_ERI3 "libint2.hpp" XTP_LIBINT_HAS_ERI3_DERIVS)
unset(CMAKE_REQUIRED_LIBRARIES)
unset(CMAKE_REQUIRED_FLAGS)

if(XTP_LIBINT_HAS_ONEBODY_DERIVS AND XTP_LIBINT_HAS_ERI2_DERIVS AND XTP_LIBINT_HAS_ERI3_DERIVS)
  message(STATUS "libint2 derivative integrals: full support found "
                 "(one-body, two-center Coulomb metric, three-center RI) "
                 "-- DFTEngine's compute_forces option will be available.")
else()
  message(WARNING
    "libint2 derivative integrals: PARTIAL OR NO support found -- "
    "one-body (overlap/kinetic/nuclear attraction): "
    "${XTP_LIBINT_HAS_ONEBODY_DERIVS}, two-center Coulomb metric: "
    "${XTP_LIBINT_HAS_ERI2_DERIVS}, three-center RI: "
    "${XTP_LIBINT_HAS_ERI3_DERIVS}. DFTEngine's compute_forces option "
    "will be unavailable (or only partially available) in this build. "
    "These are three independent libint2 configure-time flags "
    "(--enable-1body / --enable-eri2 / --enable-eri3); this is NOT a "
    "VOTCA-XTP defect, it reflects what this specific libint2 "
    "installation was actually built with.")
endif()

find_package(ecpint REQUIRED)
set_package_properties(ecpint PROPERTIES TYPE REQUIRED PURPOSE "Calculates Gaussian integrals over pseudo potentials")
message(STATUS "Found libecpint: ${ecpint_DIR}")

find_package(HDF5 REQUIRED COMPONENTS CXX)
set_package_properties(HDF5 PROPERTIES TYPE REQUIRED PURPOSE "Used to read/write HDF5 data files")

# https://github.com/votca/xtp/issues/436, hdf5-1.10.4 generates a implicitly-declared operator warning
if(HDF5_VERSION VERSION_GREATER_EQUAL 1.10.4 AND HDF5_VERSION VERSION_LESS_EQUAL 1.10.6)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-Wno-deprecated-copy" COMPILER_SUPPORTS_WNO_DEPRECATED_COPY)
  if(COMPILER_SUPPORTS_WNO_DEPRECATED_COPY)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-copy")
  endif()
endif()

option(USE_CUDA "Use a Nvidia GPU to speed up computations" OFF)
add_feature_info(USE_CUDA USE_CUDA "Use CUDA for GPU acceleration")
########################################################################
#Find external packages
########################################################################
if(USE_CUDA)
  find_package(CUDA REQUIRED)
  set_package_properties(CUDA PROPERTIES TYPE REQUIRED PURPOSE "Used for GPU acceleration")
endif(USE_CUDA)

########################################################################
# Checks what linear algebra packages are installed                    #
########################################################################



########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFileCXX)
foreach(HEADER algorithm)
  check_include_file_cxx(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

set(XTP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE STRING "XTP Source Path")

######################################
# Include the following subdirectory #
######################################
add_subdirectory(include/votca/xtp)
add_subdirectory(scripts)
add_subdirectory(src)
add_subdirectory(share)
