include(ExternalProject)
# NOTE: --enable-1body=1, --enable-eri2=1, --enable-eri3=1 were added on
# top of the original (bare --enable-eri2 --enable-eri3, no 1body flag at
# all) configuration. Without them, libint2 defaults to derivative order 0
# for all three integral types (confirmed directly against libint2 2.7.1's
# own configure.ac: the fallback branch for each of these three flags sets
# derivative order to 0/disabled when unspecified). Requesting
# deriv_order=1 on an Engine built against a libint2 compiled this way
# segfaults inside libint2's own internals (observed directly: a
# diagnostic print placed immediately after engine.compute() in
# libxtp/libint2_derivative_calls.cc never printed, meaning the crash
# happens inside compute() itself, before returning).
# =1 enables first derivatives (nuclear gradients); would need =2 if
# second derivatives (Hessians) are ever wanted later. eri2/eri3 (not just
# 1body) are included here because XTP's SCF uses RI (density-fitting,
# two- and three-center integrals), not full four-center ERIs -- the
# planned RI-J/RI-K gradient work needs derivatives of exactly these two
# integral types.
ExternalProject_Add(Libint_build
  URL https://github.com/evaleev/libint/archive/v2.7.1.tar.gz
  URL_MD5 78507c1d3de77bcd519dfd8d6afc8468
  BUILD_IN_SOURCE 1
  PREFIX libint INSTALL_DIR libint/install
  CONFIGURE_COMMAND ./autogen.sh COMMAND ./configure --prefix=<INSTALL_DIR> --libdir=<INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR} --disable-shared --with-pic --enable-eri2=1 --enable-eri3=1 --enable-1body=1 --with-max-am=6
  BUILD_BYPRODUCTS <INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR}/libint2.a
)
ExternalProject_get_property(Libint_build INSTALL_DIR)

# create emtpy directory for imported targeti below
file(MAKE_DIRECTORY "${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")
add_library(Libint2::int2 UNKNOWN IMPORTED)
set_target_properties(Libint2::int2 PROPERTIES
  IMPORTED_LOCATION "${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libint2.a"
  INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")
add_dependencies(Libint2::int2 Libint_build)

# no need to install as static libint.a get pulled into our shared lib
# install(DIRECTORY ${INSTALL_DIR}/ DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS) 
