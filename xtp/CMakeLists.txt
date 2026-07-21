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

# libint2 derivative-integral support (needed for DFTEngine's
# compute_forces option) depends on flags passed at libint2's OWN
# configure time (--enable-1body / --enable-eri2 / --enable-eri3),
# independently of each other and of whatever this libint2 install was
# otherwise built with. This was previously an actual configure-time
# check (via check_cxx_symbol_exists against libint2's own
# LIBINT_INCLUDE_ONEBODY/_ERI2/_ERI3 macros), but every attempt at that
# -- first passing imported targets directly via
# CMAKE_REQUIRED_LIBRARIES, then extracting their properties manually
# -- gave a wrong answer (false "not found") on at least one real,
# known-good machine, and one attempt caused a hard, build-breaking
# configure error on a CI runner besides. try_compile's own isolated
# sub-project, which check_cxx_symbol_exists uses internally, does not
# behave consistently enough across different libint2 discovery
# mechanisms (CONFIG-mode vs classic MODULE-mode Find<Pkg>.cmake
# scripts) and package managers to be worth the risk for what is
# ultimately a convenience message, not something anything else in the
# build depends on. A plain, static note is more honest than a check
# that has repeatedly gotten the wrong answer while looking
# authoritative. The actual, reliable answer is determined at COMPILE
# time by the per-category #ifdef guards in
# xtp/src/libxtp/libint2_derivative_calls.cc (checking these same
# macros directly, with no try_compile involved at all -- compiled as
# part of the real votca_xtp target, which has shown neither failure
# mode), and reported to the user at RUNTIME by
# DFTEngine::Initialize's own check, which fires immediately if
# compute_forces=true is requested on a build that cannot actually do
# it.
message(STATUS
  "Note: DFTEngine's compute_forces option needs this libint2 to have "
  "been configured with --enable-1body, --enable-eri2, and "
  "--enable-eri3 (independently of each other, and of whatever else "
  "this libint2 install was built with). This is not checked at "
  "configure time; if compute_forces does not work at runtime, verify "
  "these flags on the libint2 installation actually being used.")

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
