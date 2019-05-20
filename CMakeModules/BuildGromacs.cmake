include(ExternalProject)
include(GNUInstallDirs)
# ftp://ftp.gromacs.org/pub/gromacs/gromacs-XXX.tar.gz is too unstable for CI, so use Gentoo Mirror
ExternalProject_Add(Gromacs_build
  URL http://distfiles.gentoo.org/distfiles/gromacs-2018.6.tar.gz
  URL_MD5 997f2b8108197780d61b23b5eb05644a
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
)
ExternalProject_get_property(Gromacs_build INSTALL_DIR)

# create emtpy directory for imported targeti below
file(MAKE_DIRECTORY "${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")
add_library(GMX::libgromacs UNKNOWN IMPORTED)
set_target_properties(GMX::libgromacs PROPERTIES
  IMPORTED_LOCATION "${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libgromacs.so"
  INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")
add_dependencies(GMX::libgromacs Gromacs_build)
set_target_properties(GMX::libgromacs PROPERTIES GROMACS_PATH "${INSTALL_DIR}/${CMAKE_INSTALL_BINDIR}")

install(DIRECTORY ${INSTALL_DIR}/ DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS) 
