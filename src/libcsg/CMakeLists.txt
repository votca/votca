option(NOGSL "Do not use GSL, leads to reduced functionality!" OFF)
if (NOGSL)
  set(NOGSL TRUE)
else(NOGSL)
  find_package(GSL)
  if(NOT GSL_FOUND)
    message(FATAL_ERROR "GSL not found,  make sure you have also installed the gsl and it's dev package")
  endif(NOT GSL_FOUND)
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
add_library(votca_csg hgversion.h ${VOTCA_SOURCES})
target_link_libraries(votca_csg ${BOOST_LIBRARIES} ${GSL_LIBRARIES}
  ${FFTW3_LIBRARIES} ${EXPAT_LIBRARIES} ${THREAD_LIBRARIES} ${MATH_LIBRARIES})
install(TARGETS votca_csg LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_csg.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_csg.pc @ONLY)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/libvotca_csg.pc DESTINATION ${LIB}/pkgconfig)

add_executable(linktest linktest.cc)
target_link_libraries(linktest votca_csg)
