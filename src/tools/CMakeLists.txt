set(TOOLS_RST_FILES)
foreach(PROG votca_property)
  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  target_link_libraries(${PROG} votca_tools)
  #not realy useful, re-enable whenever it is ready to replace csg_property
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

  if(ENABLE_TESTING)
    add_test(integration_${PROG}Help ${PROG} --help)
    set_tests_properties(integration_${PROG}Help PROPERTIES LABELS "tools;votca;integration")
  endif(ENABLE_TESTING)

  if(VOTCA_SPHINX_DIR)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/tools/${PROG}.rst COMMAND $<TARGET_FILE:VOTCA::votca_help2doc> -n ${CMAKE_CURRENT_BINARY_DIR}/${PROG} -o ${VOTCA_SPHINX_DIR}/tools/${PROG}.rst
      DEPENDS $<TARGET_FILE:VOTCA::votca_help2doc> ${PROG})
    list(APPEND TOOLS_RST_FILES ${VOTCA_SPHINX_DIR}/tools/${PROG}.rst)
  endif(VOTCA_SPHINX_DIR)
endforeach(PROG)

if(VOTCA_SPHINX_DIR)
  add_custom_target(doc-tools-prog DEPENDS ${TOOLS_RST_FILES})
  add_dependencies(doc-tools doc-tools-prog)
endif()
