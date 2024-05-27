set(XTP_RST_FILES)
foreach(PROG xtp_map xtp_run xtp_tools xtp_parallel)
  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  target_link_libraries(${PROG} votca_xtp)
  install(TARGETS ${PROG} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  if (BUILD_MANPAGES)
    add_custom_command(OUTPUT ${PROG}.man
      COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_help2doc> --name ${PROG} --format groff --out ${PROG}.man
      COMMENT "Building ${PROG} manpage"
      DEPENDS VOTCA::votca_help2doc ${PROG})
    add_custom_target(${PROG}_manpage DEPENDS ${PROG}.man)
    add_dependencies(manpages ${PROG}_manpage)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROG}.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME ${PROG}.1)
  endif (BUILD_MANPAGES)

  if(SPHINX_FOUND)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/xtp/${PROG}.rst
      COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_help2doc> --name ${PROG} --format rst --out ${VOTCA_SPHINX_DIR}/xtp/${PROG}.rst
      COMMENT "Building ${PROG} rst doc"
      DEPENDS VOTCA::votca_help2doc ${PROG})
    list(APPEND XTP_RST_FILES ${VOTCA_SPHINX_DIR}/xtp/${PROG}.rst)
  endif()

  add_test(integration_${PROG}Help ${PROG} --help)
endforeach(PROG)

if(SPHINX_FOUND)
  add_custom_target(xtp-tools-rst DEPENDS ${XTP_RST_FILES})
  add_dependencies(sphinx-deps xtp-tools-rst)
endif()
