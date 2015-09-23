configure_file(help2t2t.in ${CMAKE_CURRENT_BINARY_DIR}/help2t2t.out @ONLY)
add_custom_target(help2t2t_build DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/help2t2t)
add_custom_command(OUTPUT help2t2t ALL COMMAND ${CMAKE_COMMAND}
    -DINPUT="help2t2t.out" -DOUTPUT="help2t2t"
    -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
    -DMERCURIAL_EXECUTABLE="${MERCURIAL_EXECUTABLE}"
    -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${CMAKE_MODULE_PATH}/gitscript.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES help2t2t)

foreach(SCRIPT csg_call csg_inverse multi_g_rdf)
  configure_file(${SCRIPT}.in ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.out @ONLY)
  add_custom_target(${SCRIPT}_build ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT})
  add_custom_command(OUTPUT ${SCRIPT} COMMAND ${CMAKE_COMMAND}
    -DINPUT="${SCRIPT}.out" -DOUTPUT="${SCRIPT}"
    -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
    -DMERCURIAL_EXECUTABLE="${MERCURIAL_EXECUTABLE}"
    -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${CMAKE_MODULE_PATH}/gitscript.cmake)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${SCRIPT})
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} DESTINATION bin)
  if (TXT2TAGS_FOUND AND BASH)
    add_custom_command(OUTPUT ${SCRIPT}.man 
      COMMAND VOTCASHARE=${CMAKE_SOURCE_DIR}/share ${BASH} ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} --help > ${SCRIPT}.help
      COMMAND ${BASH} ${CMAKE_CURRENT_BINARY_DIR}/help2t2t ${SCRIPT}.help > ${SCRIPT}.t2t 
      COMMAND ${TXT2TAGS_EXECUTABLE} -q -t man -i ${SCRIPT}.t2t -o ${SCRIPT}.man
      DEPENDS help2t2t_build ${SCRIPT})
    add_custom_target(${SCRIPT}_manpage DEPENDS ${SCRIPT}.man)
    add_dependencies(manpages ${SCRIPT}_manpage)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.man DESTINATION ${MAN}/man1 RENAME ${SCRIPT}.1)
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${SCRIPT}.help ${SCRIPT}.t2t)
  endif(TXT2TAGS_FOUND AND BASH)
endforeach(SCRIPT)

option(WITH_RC_FILES "Install votca rc files, no need when installing under /usr" ON)
install(FILES csg-completion.bash DESTINATION ${DATA}/rc)
