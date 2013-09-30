find_package(Doxygen)
if (DOXYGEN_FOUND)
  configure_file(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(html COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile 
    COMMENT "Build doxygen documentation")
endif (DOXYGEN_FOUND)
