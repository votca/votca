include(ExternalProject)
ExternalProject_Add(Libint_build
  URL https://github.com/evaleev/libint/archive/v2.7.1.tar.gz
  URL_MD5 78507c1d3de77bcd519dfd8d6afc8468
  BUILD_IN_SOURCE 1
  PREFIX libint INSTALL_DIR libint/install
  CONFIGURE_COMMAND ./autogen.sh COMMAND ./configure --prefix=<INSTALL_DIR> --libdir=<INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR} --disable-shared --with-pic --enable-eri2 --enable-eri3 --with-max-am=6
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
