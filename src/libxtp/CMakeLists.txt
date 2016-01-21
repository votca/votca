option(WITH_GSL "Use GSL, disabling leads to reduced functionality!" OFF)
if (WITH_GSL)
  find_package(GSL)
  if(NOT GSL_FOUND)
    message(FATAL_ERROR "GSL not found,  make sure you have also installed the gsl and it's dev package")
  endif(NOT GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  set(GSL_PKG "gsl")
else(WITH_GSL)
  set(NOGSL TRUE)
  set(GSL_PKG)
endif(WITH_GSL)

option(WITH_MKL "Checking if MKL is installed in order to improve UBLAS performance" ON)
if (WITH_MKL)
  find_package(MKL)
  if(NOT MKL_FOUND AND CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    #TODO: try to append -mkl to cflags
  endif(NOT MKL_FOUND AND CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
endif(WITH_MKL)

configure_file(votca_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(hgversion-xtp COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -DMERCURIAL_EXECUTABLE="${MERCURIAL_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/hgversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES hgversion.h)


file(GLOB_RECURSE VOTCA_SOURCES *.cc)
file(GLOB_RECURSE NOT_VOTCA_SOURCES version_nb.cc segmenttype.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})

add_library(votca_xtp  ${VOTCA_SOURCES})
add_dependencies(votca_xtp hgversion-xtp)
target_link_libraries(votca_xtp ${VOTCA_CSG_LIBRARIES} ${MKL_LIBRARIES} ${GSL_LIBRARIES}  votca_xtp_moo ${VOTCA_TOOLS_LIBRARIES} ${BOOST_LIBRARIES} ${LIBXC_LIBRARIES})
install(TARGETS votca_xtp LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_xtp.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_xtp.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libvotca_xtp.pc DESTINATION ${LIB}/pkgconfig)
