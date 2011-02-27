find_package(TXT2TAGS)
if (TXT2TAGS_FOUND)
  add_custom_command(OUTPUT votca-tools.man 
    COMMAND ${TXT2TAGS_EXECUTABLE} -t man -i ${CMAKE_CURRENT_SOURCE_DIR}/votca-tools.t2t -o votca-tools.man
    COMMENT "Building votca_tools manpage")
  add_custom_target(manpages DEPENDS votca-tools.man)
  add_dependencies(votca_tools manpages)
  install(FILES votca_tools.man DESTINATION ${MAN}/man7)
endif (TXT2TAGS_FOUND)

