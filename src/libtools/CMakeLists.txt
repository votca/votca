find_package(EXPAT REQUIRED)
include_directories(${EXPAT_INCLUDE_DIRS})

option(NOFFTW "Do not use FFTW3" OFF)
if (NOFFTW)
  set(NOFFTW TRUE)
else(NOFFTW)
  find_package(FFTW3 REQUIRED)
  include_directories(${FFTW3_INCLUDE_DIR})
endif(NOFFTW)

option(NOGSL "Do not use GSL" OFF)
if (NOGSL)
  set(NOGSL TRUE)
else(NOGSL)
  find_package(GSL REQUIRED)
  include_directories(${GSL_INCLUDE_DIR})
endif(NOGSL)

add_definitions( -DHAVE_CONFIG_H )
configure_file(config.h.in config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

if (MERCURIAL_FOUND)
  MERCURIAL_HG_INFO(${PROJECT_SOURCE_DIR} TOOLS)
  MESSAGE("Current revision is ${TOOLS_HG_ID}")
else(MERCURIAL_FOUND)
  set (TOOLS_HG_ID)
endif(MERCURIAL_FOUND)
set (HG_HEADER "hgversion.h")
add_custom_command(OUTPUT ${HG_HEADER}
  COMMAND ${CMAKE_COMMAND} -E echo 
  static const std::string hgversion = "\\\"hgid: ${TOOLS_HG_ID}\\\"\;" > ${HG_HEADER})

file(GLOB_RECURSE VOTCA_SOURCES *.cc)
file(GLOB_RECURSE NOT_VOTCA_SOURCES version_nb.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})
add_library(votca_tools hgversion.h ${VOTCA_SOURCES})
target_link_libraries(votca_tools ${BOOST} 
  ${FFTW3_LIBRARIES} ${GSL_LIBRARIES} ${EXPAT_LIBRARIES})
install(TARGETS votca_tools LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_tools.pc.in libvotca_tools.pc @ONLY)
install(FILES libvotca_tools.pc DESTINATION ${LIB}/pkgconfig)
