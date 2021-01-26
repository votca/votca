foreach(SCRIPT xtp_update_mapfile xtp_basisset xtp_makeauxbasis xtp_autogen_mapping xtp_qmmm2qm)
  configure_file(${SCRIPT}.in ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.out @ONLY)
  add_custom_target(${SCRIPT}_build ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT})
  add_custom_command(OUTPUT ${SCRIPT} COMMAND ${CMAKE_COMMAND}
    -DINPUT="${SCRIPT}.out" -DOUTPUT="${SCRIPT}"
    -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
    -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${PROJECT_SOURCE_DIR}/CMakeModules/gitscript.cmake)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${SCRIPT})
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach(SCRIPT)
