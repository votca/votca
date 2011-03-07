file(GLOB VOTCA_SOURCES *.cc *.cpp)
#file(GLOB NOT_VOTCA_SOURCES xxx.cc)
#list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})
add_library(votca_kmc  ${VOTCA_SOURCES})
target_link_libraries(votca_kmc ${VOTCA_TOOLS_LIBRARIES})
install(TARGETS votca_kmc LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_kmc.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_kmc.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libvotca_kmc.pc DESTINATION ${LIB}/pkgconfig)
