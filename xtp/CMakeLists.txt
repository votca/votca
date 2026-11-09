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

# libint2 derivative integrals (needed for the compute_forces option of the
# DFT engine) require libint2 configured with --enable-1body, --enable-eri2
# and --enable-eri3. Checked by reading libint2/config.h;
# libint2_derivative_calls.cc checks the same macros at compile time.
include(CheckLibintFeatures)
votca_check_libint_features()
add_feature_info(LIBINT_FORCES VOTCA_LIBINT_FORCES "libint2 derivative integrals for analytic DFT forces")

find_package(ecpint REQUIRED)
set_package_properties(ecpint PROPERTIES TYPE REQUIRED PURPOSE "Calculates Gaussian integrals over pseudo potentials")
message(STATUS "Found libecpint: ${ecpint_DIR}")

find_package(HDF5 REQUIRED COMPONENTS CXX)
set_package_properties(HDF5 PROPERTIES TYPE REQUIRED PURPOSE "Used to read/write HDF5 data files")

# FragmentSaturator relaxes newly placed H-saturation atoms with
# OpenBabel (MMFF94/UFF). Required rather than optional, so that results do
# not depend on which optional packages a build has.
#
# Known upstream bug (openbabel/openbabel #2305, #2374, #2264):
# OpenBabel3Config.cmake can set OpenBabel3_FOUND to FALSE because it
# references the imported targets maeparser and coordgen, which many installs
# lack (VOTCA needs neither). If find_package fails, look for "targets are
# referenced, but are missing: maeparser coordgen" in the configure output.
find_package(OpenBabel3 REQUIRED)
set_package_properties(OpenBabel3 PROPERTIES TYPE REQUIRED PURPOSE "Places and relaxes new H atoms when saturating covalently-bonded fragment boundaries")

# https://github.com/votca/xtp/issues/436, hdf5-1.10.4 generates a implicitly-declared operator warning
if(HDF5_VERSION VERSION_GREATER_EQUAL 1.10.4 AND HDF5_VERSION VERSION_LESS_EQUAL 1.10.6)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-Wno-deprecated-copy" COMPILER_SUPPORTS_WNO_DEPRECATED_COPY)
  if(COMPILER_SUPPORTS_WNO_DEPRECATED_COPY)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-copy")
  endif()
endif()

if(USE_CUDA)
  message(FATAL_ERROR "USE_CUDA: GPU offloading with CUDA has been removed "
    "from xtp; configure without -DUSE_CUDA=ON.")
endif()

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
