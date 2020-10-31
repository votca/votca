include(ExternalProject)
include(GNUInstallDirs)
enable_language(C)
# ftp://ftp.gromacs.org/pub/gromacs/gromacs-XXX.tar.gz is too unstable for CI, so use Gentoo Mirror
ExternalProject_Add(Gromacs_build
  URL http://distfiles.gentoo.org/distfiles/gromacs-2019.6.tar.gz https://ftp.gromacs.org/pub/gromacs/gromacs-2019.6.tar.gz
  URL_MD5 a9c41f5b9b11f738a6e04c2427206e09
  PREFIX gromacs INSTALL_DIR gromacs/install
  CMAKE_ARGS 
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DGMX_MPI:BOOL=OFF
    -DGMX_THREAD_MPI:BOOL=ON 
    -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
  BUILD_BYPRODUCTS <INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR}/libgromacs.so
)
ExternalProject_get_property(Gromacs_build INSTALL_DIR)

# create emtpy directory for imported targeti below
file(MAKE_DIRECTORY "${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")
add_library(Gromacs::libgromacs UNKNOWN IMPORTED)
set_target_properties(Gromacs::libgromacs PROPERTIES
  IMPORTED_LOCATION "${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libgromacs.so"
  INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR}"
  INTERFACE_COMPILE_OPTIONS "-DGMX_DOUBLE=0")
add_dependencies(Gromacs::libgromacs Gromacs_build)
set_target_properties(Gromacs::libgromacs PROPERTIES GROMACS_PATH "${INSTALL_DIR}/${CMAKE_INSTALL_BINDIR}")

install(DIRECTORY ${INSTALL_DIR}/ DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS) 
