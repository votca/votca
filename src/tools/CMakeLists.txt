set(XTP_PROGS) 
foreach(PROG xtp_map xtp_run xtp_tools xtp_parallel xtp_dump)
  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)

  if( ${PROG} STREQUAL "xtp_map" )
     add_executable(xtp_map xtp_map.cc Md2QmEngine.cc)
  else()
     add_executable(${PROG} ${${PROG}_SOURCES})
  endif()

  target_link_libraries(${PROG} votca_xtp)
  install(TARGETS ${PROG} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  if (BUILD_MANPAGES)
     add_custom_command(OUTPUT ${PROG}.man
       COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PROG} --man > ${PROG}.man
       COMMENT "Building ${PROG} manpage"
       DEPENDS ${PROG})
     add_custom_target(${PROG}_manpage DEPENDS ${PROG}.man)
     add_dependencies(manpages ${PROG}_manpage)
     install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROG}.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME ${PROG}.1)
  endif (BUILD_MANPAGES)
  list(APPEND XTP_PROGS "${PROG}")

  if(ENABLE_TESTING)
    add_test(${PROG}Help ${PROG} --help)
    # run tests for tools, csg, moo, ctp as well for coverage
    set_tests_properties(${PROG}Help PROPERTIES LABELS "xtp;tools;csg;ctp;moo;votca")
  endif(ENABLE_TESTING)
endforeach(PROG)
set(XTP_PROGS "${XTP_PROGS}" PARENT_SCOPE)
