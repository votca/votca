include(ExternalProject)
include(GNUInstallDirs)
ExternalProject_Add(Libint_build
  URL https://github.com/evaleev/libint/archive/v2.6.0.tar.gz
  URL_MD5 7241fc0005f37394302083c55a36e0e6
  BUILD_IN_SOURCE 1
  PREFIX libint INSTALL_DIR libint/install
  CONFIGURE_COMMAND ./autogen.sh COMMAND ./configure --prefix=<INSTALL_DIR> --libdir=<INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR} --disable-shared --with-pic --enable-eri2 --enable-eri3 --with-max-am=6
  BUILD_BYPRODUCTS <INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR}/libint2.a
)
ExternalProject_get_property(Libint_build INSTALL_DIR)

# create emtpy directory for imported targeti below
file(MAKE_DIRECTORY "${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")
add_library(PkgConfig::LIBINT UNKNOWN IMPORTED)
set_target_properties(PkgConfig::LIBINT PROPERTIES
  IMPORTED_LOCATION "${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libint2.a"
  INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")
add_dependencies(PkgConfig::LIBINT Libint_build)

# no need to install as static libint.a get pulled into our shared lib
# install(DIRECTORY ${INSTALL_DIR}/ DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS) 
