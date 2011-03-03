configure_file(votca_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

file(GLOB_RECURSE VOTCA_SOURCES *.cc)
file(GLOB NOT_VOTCA_SOURCES version_nb.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})
add_library(votca_md2qm  ${VOTCA_SOURCES})
target_link_libraries(votca_md2qm ${VOTCA_CSG_LIBRARIES} ${VOTCA_KMC_LIBRARIES} ${VOTCA_MOO_LIBRARIES} ${VOTCA_TOOLS_LIBRARIES})
install(TARGETS votca_md2qm LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_md2qm.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_md2qm.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libvotca_md2qm.pc DESTINATION ${LIB}/pkgconfig)
