foreach(PROG template_serial template_threaded)
  file(GLOB ${PROG}_SOURCES ${PROG}.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  target_link_libraries(${PROG} votca_csg)
  install(FILES ${${PROG}_SOURCES} DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/template)
 
  if(ENABLE_TESTING)
    add_test(${PROG}Help ${PROG} --help)
    # run tests for tools (for coverage) as well
    set_tests_properties(${PROG}Help PROPERTIES LABELS "csg;tools;votca")
  endif(ENABLE_TESTING)
endforeach(PROG)

install(FILES CMakeLists.txt.template DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/template RENAME CMakeLists.txt)
