option(WITH_GSL "Use GSL, disabling leads to reduced functionality!" ON)
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

configure_file(votca_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(hgversion COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${CMAKE_MODULE_PATH}/hgversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES hgversion.h)


file(GLOB_RECURSE VOTCA_SOURCES *.cc)
file(GLOB_RECURSE NOT_VOTCA_SOURCES version_nb.cc segmenttype.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})

add_library(votca_ctp  ${VOTCA_SOURCES})
add_dependencies(votca_ctp hgversion)
target_link_libraries(votca_ctp ${VOTCA_CSG_LIBRARIES}  ${GSL_LIBRARIES} ${VOTCA_KMCOLD_LIBRARIES} ${VOTCA_MOO_LIBRARIES} ${VOTCA_TOOLS_LIBRARIES} ${BOOST_LIBRARIES})
install(TARGETS votca_ctp LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_ctp.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_ctp.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libvotca_ctp.pc DESTINATION ${LIB}/pkgconfig)
