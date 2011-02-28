if (TXT2TAGS_FOUND)
  add_custom_command(OUTPUT votca-csg.man 
    COMMAND ${TXT2TAGS_EXECUTABLE} -q -t man -i ${CMAKE_CURRENT_SOURCE_DIR}/votca-csg.t2t -o votca-csg.man
    COMMENT "Building votca-csg manpage")
  add_custom_target(manpages DEPENDS votca-csg.man)
  add_dependencies(votca_csg manpages)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/votca-csg.man DESTINATION ${MAN}/man7)
endif (TXT2TAGS_FOUND)

