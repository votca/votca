configure_file(votca_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(hgversion COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${CMAKE_MODULE_PATH}/hgversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES hgversion.h)


file(GLOB_RECURSE VOTCA_SOURCES *.cc)
file(GLOB NOT_VOTCA_SOURCES version_nb.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})
add_library(votca_md2qm  ${VOTCA_SOURCES})
add_dependencies(votca_md2qm hgversion)
target_link_libraries(votca_md2qm ${VOTCA_CSG_LIBRARIES} ${VOTCA_KMCOLD_LIBRARIES} ${VOTCA_MOO_LIBRARIES} ${VOTCA_TOOLS_LIBRARIES} ${SQLITE3_LIBRARIES})
install(TARGETS votca_md2qm LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_md2qm.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_md2qm.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libvotca_md2qm.pc DESTINATION ${LIB}/pkgconfig)
