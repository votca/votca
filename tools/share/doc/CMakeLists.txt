if(VOTCA_SPHINX_DIR)
  file(GLOB RSTS ${CMAKE_CURRENT_SOURCE_DIR}/*.rst)
  set(MANUAL_RST_FILES)
  foreach(_RST ${RSTS})
    get_filename_component(_FILE "${_RST}" NAME)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/tools/${_FILE}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_RST} ${VOTCA_SPHINX_DIR}/tools/${_FILE}
      DEPENDS ${_RST})
    list(APPEND MANUAL_RST_FILES ${VOTCA_SPHINX_DIR}/tools/${_FILE})
  endforeach()
  add_custom_target(doc-tools-manual DEPENDS ${MANUAL_RST_FILES})
  add_dependencies(doc-tools doc-tools-manual)
endif()
