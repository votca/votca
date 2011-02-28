find_package(EXPAT REQUIRED)
include_directories(${EXPAT_INCLUDE_DIR})
set (EXPAT_CFLAGS_PKG "-I${EXPAT_INCLUDE_DIRS}")
string(REGEX REPLACE "^(.*)/lib([^/]*)\\.[^.]*$" "-L\\1 -l\\2" EXPAT_LIBS_PKG ${EXPAT_LIBRARY})

find_package(PkgConfig)
option(NOFFTW "Do not use FFTW3" OFF)
if (NOFFTW)
  set(NOFFTW TRUE)
else(NOFFTW)
  find_package(FFTW3 REQUIRED)
  include_directories(${FFTW3_INCLUDE_DIRS})
endif(NOFFTW)

option(NOGSL "Do not use GSL" OFF)
if (NOGSL)
  set(NOGSL TRUE)
else(NOGSL)
  find_package(GSL REQUIRED)
  include_directories(${GSL_INCLUDE_DIRS})
endif(NOGSL)

configure_file(votca_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)
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
file(GLOB_RECURSE NOT_VOTCA_SOURCES version_nb.cc test.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})
add_library(votca_tools hgversion.h ${VOTCA_SOURCES})
target_link_libraries(votca_tools ${BOOST_LIBRARIES} ${GSL_LIBRARIES}
  ${FFTW3_LIBRARIES} ${EXPAT_LIBRARIES} ${THREAD_LIBRARIES} ${MATH_LIBRARIES})
install(TARGETS votca_tools LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_tools.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_tools.pc @ONLY)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/libvotca_tools.pc DESTINATION ${LIB}/pkgconfig)

add_executable(test test.cc)
target_link_libraries(test votca_tools)
