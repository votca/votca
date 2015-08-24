include_directories(${CMAKE_SOURCE_DIR}/include)

foreach(PROG template_serial template_threaded)
  file(GLOB ${PROG}_SOURCES ${PROG}.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  target_link_libraries(${PROG} votca_csg)
 install(FILES ${${PROG}_SOURCES} DESTINATION ${DATA}/template)
endforeach(PROG)

install(FILES CMakeLists.txt.template DESTINATION ${DATA}/template RENAME CMakeLists.txt)
install(FILES Makefile.template DESTINATION ${DATA}/template RENAME Makefile)
