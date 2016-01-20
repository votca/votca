find_package(Doxygen)
if (DOXYGEN_FOUND)
  configure_file(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    COMMENT "Build doxygen documentation")
  add_custom_target(html-xtp DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
  if(NOT TARGET html)
    add_custom_target(html)
  endif()
  add_dependencies(html html-xtp)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES html)
endif (DOXYGEN_FOUND)
