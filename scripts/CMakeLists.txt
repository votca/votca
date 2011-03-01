configure_file(help2t2t.in ${CMAKE_CURRENT_BINARY_DIR}/help2t2t @ONLY)

foreach(SCRIPT csg_call csg_inverse multi_g_rdf)
  configure_file(${SCRIPT}.in ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} @ONLY)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} DESTINATION bin)
  if (TXT2TAGS_FOUND AND BASH)
    add_custom_command(OUTPUT ${SCRIPT}.man 
      COMMAND CSGSHARE=${CMAKE_SOURCE_DIR}/share ${BASH} ${SCRIPT} --help > ${SCRIPT}.help 
      COMMAND ${BASH} ${CMAKE_CURRENT_BINARY_DIR}/help2t2t ${SCRIPT}.help > ${SCRIPT}.t2t 
      COMMAND ${TXT2TAGS_EXECUTABLE} -q -t man -i ${SCRIPT}.t2t -o ${SCRIPT}.man
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} ${CMAKE_CURRENT_BINARY_DIR}/help2t2t)
    add_custom_target(${SCRIPT}_manpage DEPENDS ${SCRIPT}.man)
    add_dependencies(manpages ${SCRIPT}_manpage)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.man DESTINATION ${MAN}/man1)
  endif(TXT2TAGS_FOUND AND BASH)
endforeach(SCRIPT)

foreach(RC csg.rc.bash csg.rc.csh)
  configure_file(${RC}.in ${CMAKE_CURRENT_BINARY_DIR}/${RC} @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${RC} DESTINATION @DATA@/rc)
endforeach(RC)

install(FILES csg-comletion.bash DESTINATION @DATA@/rc)
