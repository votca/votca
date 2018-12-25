if (TXT2TAGS_FOUND)
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/votca-csg.man 
    COMMAND ${TXT2TAGS_EXECUTABLE} -q -t man -i ${CMAKE_CURRENT_SOURCE_DIR}/votca-csg.t2t -o ${CMAKE_CURRENT_BINARY_DIR}/votca-csg.man
    COMMENT "Building votca-csg manpage")
  add_custom_target(votca-csg_manpage DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/votca-csg.man)
  add_dependencies(manpages votca-csg_manpage)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/votca-csg.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man7 RENAME votca-csg.7)
endif (TXT2TAGS_FOUND)

