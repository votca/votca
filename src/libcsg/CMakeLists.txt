option(NOGMX "Do not build gromacs reader/writer, leads to reduced functionality!" OFF)

file(GLOB GMX_SOURCES modules/io/gmx*.cc modules/io/gro*.cc modules/io/pdbtopology*.cc)
file(GLOB NOT_GMX_SOURCES modules/io/gmx_print_version.cc modules/io/gmx_version_nb.cc)
list(REMOVE_ITEM GMX_SOURCES ${NOT_GMX_SOURCES})

file(GLOB IO_SOURCES modules/io/*.cc)
file(GLOB NOT_IO_SOURCES modules/io/mdp*.cc)
list(REMOVE_ITEM IO_SOURCES ${GMX_SOURCES} ${NOT_GMX_SOURCES} ${NOT_IO_SOURCES})

if (NOGMX)
  set(GMX)
  set(GMX_SOURCES)
else(NOGMX)
  #prefer libgromacs over libgmx
  find_package(GROMACS)
  if (GROMACS_FOUND)
    set(GMX_LIBS "${GROMACS_LIBRARIES}")
    include_directories(${GROMACS_INCLUDE_DIRS})
    set(GMX 50)
  else(GROMACS_FOUND)
    find_package(GMX)
    if(GMX_FOUND)
      set(GMX_LIBS "${GMX_LIBRARIES}")
      include_directories(${GMX_INCLUDE_DIRS})
      set(GMX 45)
    else(GMX_FOUND)
      message(FATAL_ERROR "gromacs not found,  make sure you have also installed the gromacs and it's dev package or disable it with -DNOGMX=ON")
    endif(GMX_FOUND)
  endif(GROMACS_FOUND)
  add_executable(gmx_print_version modules/io/gmx_print_version.cc)
  target_link_libraries(gmx_print_version ${GMX_LIBS})
  add_custom_command(OUTPUT gmx_libs_version.h COMMAND gmx_print_version > gmx_libs_version.h DEPENDS gmx_print_version) 
  list(APPEND GMX_SOURCES gmx_libs_version.h)
endif(NOGMX)

configure_file(votca_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set (HG_HEADER "hgversion.h")
add_custom_command(OUTPUT ${HG_HEADER}
  COMMAND ${CMAKE_COMMAND} -E echo 
  static const std::string hgversion = "\\\"${CSG_HG_ID}\\\"\;" > ${HG_HEADER})

file(GLOB VOTCA_SOURCES *.cc)
file(GLOB NOT_VOTCA_SOURCES version_nb.cc test.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})
add_library(votca_csg hgversion.h ${VOTCA_SOURCES} ${GMX_SOURCES} ${IO_SOURCES})
target_link_libraries(votca_csg ${VOTCA_TOOLS_LIBRARIES} ${GMX_LIBS})
install(TARGETS votca_csg LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_csg.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_csg.pc @ONLY)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/libvotca_csg.pc DESTINATION ${LIB}/pkgconfig)

add_executable(linktest linktest.cc)
target_link_libraries(linktest votca_csg)
