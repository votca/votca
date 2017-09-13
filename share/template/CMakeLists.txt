include_directories(${CMAKE_SOURCE_DIR}/include)

foreach(PROG template_serial template_threaded)
  file(GLOB ${PROG}_SOURCES ${PROG}.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  target_link_libraries(${PROG} votca_csg)
 install(FILES ${${PROG}_SOURCES} DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/template)
endforeach(PROG)

install(FILES CMakeLists.txt.template DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/template RENAME CMakeLists.txt)
install(FILES Makefile.template DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/template RENAME Makefile)
