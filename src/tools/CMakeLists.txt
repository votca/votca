foreach(PROG votca_property)
  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  target_link_libraries(${PROG} votca_tools)
  #not realy useful, re-enable whenever it is ready to replace csg_property
  #install(TARGETS ${PROG} RUNTIME DESTINATION bin)

  #help2t2t is part of csg for now, re-enable whenever votca_tools is ready
  # or if votca issue 94 and votca-ctp issue 16 is solved
  #if (TXT2TAGS_FOUND AND BASH)
  #  add_custom_command(OUTPUT ${PROG}.man
  #    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PROG} --help > ${PROG}.help
  #    COMMAND ${BASH} ${CMAKE_BINARY_DIR}/scripts/help2t2t ${PROG}.help > ${PROG}.t2t
  #    COMMAND ${TXT2TAGS_EXECUTABLE} -q -t man -i ${PROG}.t2t -o ${PROG}.man
  #    DEPENDS scripts/help2t2t ${PROG})
  #  add_custom_target(${PROG}_manpage DEPENDS ${PROG}.man)
  #  add_dependencies(manpages ${PROG}_manpage)
  #  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROG}.man DESTINATION ${MAN}/man1 RENAME ${PROG}.1)
  #  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${PROG}.help ${PROG}.t2t)
  #endif(TXT2TAGS_FOUND AND BASH)
endforeach(PROG)
