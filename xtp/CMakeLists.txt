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

# Real, hands-on history behind this check: LIBINT2_MAX_DERIV_ORDER (a
# macro libint2's own headers define) is library-wide -- the maximum
# derivative order ANY operator supports -- but individual libint2
# builds/packages have shown, on real CI runs, that it can report
# LIBINT2_MAX_DERIV_ORDER>=1 while derivative support for SPECIFIC
# operators (nuclear attraction with a point charge; even plain overlap
# on one architecture) either silently returns zero-filled buffers or
# segfaults outright inside engine.compute() itself -- a crash that
# cannot be caught, detected, or guarded against from ANY
# application-level code once it happens. The only reliable way to catch
# this is to actually run a real deriv_order=1 computation before the
# main build exists at all, so a broken libint2 installation is
# diagnosed at configure time rather than as a segfault deep inside a
# real user's SCF run.
try_run(
  XTP_LIBINT2_DERIVATIVES_RUN_RESULT
  XTP_LIBINT2_DERIVATIVES_COMPILE_RESULT
  ${CMAKE_BINARY_DIR}/CMakeTmp/xtp_libint2_derivatives_check
  ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/check_libint2_derivatives.cc
  LINK_LIBRARIES Libint2::int2
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  RUN_OUTPUT_VARIABLE XTP_LIBINT2_DERIVATIVES_OUTPUT
)

if(XTP_LIBINT2_DERIVATIVES_COMPILE_RESULT AND
   "${XTP_LIBINT2_DERIVATIVES_RUN_RESULT}" STREQUAL "0")
  message(STATUS "libint2 derivative integrals: working "
                 "(overlap, kinetic, nuclear attraction, Coulomb metric "
                 "all verified at configure time)")
  set(XTP_LIBINT2_DERIVATIVES_WORK TRUE)
else()
  message(WARNING
    "libint2 derivative integrals: NOT usable on this build "
    "(compile result: ${XTP_LIBINT2_DERIVATIVES_COMPILE_RESULT}, "
    "run result: ${XTP_LIBINT2_DERIVATIVES_RUN_RESULT}). "
    "Output:\n${XTP_LIBINT2_DERIVATIVES_OUTPUT}\n"
    "Analytic nuclear forces (DFTEngine's compute_forces option) will "
    "be unavailable in this build -- rebuild libint2 with working "
    "derivative-integral support (--enable-1body=1 --enable-eri2=1 at "
    "libint2 configure time) to enable this feature. This is NOT a "
    "VOTCA-XTP code defect; it reflects what this specific libint2 "
    "installation actually supports at runtime, which the "
    "LIBINT2_MAX_DERIV_ORDER header macro alone was found insufficient "
    "to detect reliably.")
  set(XTP_LIBINT2_DERIVATIVES_WORK FALSE)
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
