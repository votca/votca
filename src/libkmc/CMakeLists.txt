file(GLOB VOTCA_SOURCES *.cc *.cpp)

configure_file(votca_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(hgversion COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -DMERCURIAL_EXECUTABLE="${MERCURIAL_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/hgversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES hgversion.h)

file(GLOB_RECURSE VOTCA_SOURCES *.cc)
file(GLOB_RECURSE NOT_VOTCA_SOURCES version_nb.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})

add_library(votca_xtp_kmc  ${VOTCA_SOURCES})
add_dependencies(votca_xtp_kmc hgversion)
target_link_libraries(votca_xtp_kmc ${VOTCA_TOOLS_LIBRARIES})
install(TARGETS votca_xtp_kmc LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_xtp_kmc.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_xtp_kmc.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libvotca_xtp_kmc.pc DESTINATION ${LIB}/pkgconfig)
