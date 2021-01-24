set(XTP_RST_FILES)
foreach(PROG xtp_map xtp_run xtp_tools xtp_parallel)
  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  target_link_libraries(${PROG} votca_xtp)
  install(TARGETS ${PROG} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  if (BUILD_MANPAGES)
    add_custom_command(OUTPUT ${PROG}.man
      COMMAND $<TARGET_FILE:VOTCA::votca_help2doc> --name $<TARGET_FILE:VOTCA::${PROG}> --format groff --out ${PROG}.man
      COMMENT "Building ${PROG} manpage"
      DEPENDS $<TARGET_FILE:VOTCA::votca_help2doc> ${PROG})
    add_custom_target(${PROG}_manpage DEPENDS ${PROG}.man)
    add_dependencies(manpages ${PROG}_manpage)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROG}.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME ${PROG}.1)
  endif (BUILD_MANPAGES)

  if(VOTCA_SPHINX_DIR)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/xtp/${PROG}.rst
      COMMAND $<TARGET_FILE:VOTCA::votca_help2doc> --name $<TARGET_FILE:VOTCA::${PROG}> --format rst --out ${VOTCA_SPHINX_DIR}/xtp/${PROG}.rst
      COMMENT "Building ${PROG} rst doc"
      DEPENDS $<TARGET_FILE:VOTCA::votca_help2doc> ${PROG})
    list(APPEND XTP_RST_FILES ${VOTCA_SPHINX_DIR}/xtp/${PROG}.rst)
  endif(VOTCA_SPHINX_DIR)

  if(ENABLE_TESTING)
    add_test(integration_${PROG}Help ${PROG} --help)
    # run tests for tools, csg as well for coverage
    set_tests_properties(integration_${PROG}Help PROPERTIES LABELS "xtp;tools;csg;votca;integration")
  endif(ENABLE_TESTING)
endforeach(PROG)

if(VOTCA_SPHINX_DIR)
  add_custom_target(doc-xtp-rst-progs DEPENDS ${XTP_RST_FILES})
  add_dependencies(doc-xtp doc-xtp-rst-progs)
endif()
