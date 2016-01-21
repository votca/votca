file(GLOB VOTCA_SOURCES *.cc *.cpp)
#file(GLOB NOT_VOTCA_SOURCES xxx.cc)
#list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})
add_library(votca_xtp_moo  ${VOTCA_SOURCES})
target_link_libraries(votca_xtp_moo ${VOTCA_TOOLS_LIBRARIES})
install(TARGETS votca_xtp_moo LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_xtp_moo.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_xtp_moo.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libvotca_xtp_moo.pc DESTINATION ${LIB}/pkgconfig)
