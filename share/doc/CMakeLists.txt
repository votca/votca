find_package(Doxygen)
if (DOXYGEN_FOUND)
  configure_file(Doxyfile.in Doxyfile @ONLY)
  add_custom_target(html COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile)
  install(DIRECTORY html DESTINATION ${DATA})
endif (DOXYGEN_FOUND)

