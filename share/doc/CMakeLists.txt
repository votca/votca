if(VOTCA_SPHINX_DIR)
  file(GLOB RSTS ${CMAKE_CURRENT_SOURCE_DIR}/*.rst)
  set(MANUAL_RST_FILES)
  foreach(_RST ${RSTS})
    get_filename_component(_FILE "${_RST}" NAME)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/xtp/${_FILE}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_RST} ${VOTCA_SPHINX_DIR}/xtp/${_FILE}
      DEPENDS ${_RST})
    list(APPEND MANUAL_RST_FILES ${VOTCA_SPHINX_DIR}/xtp/${_FILE})
  endforeach()

  find_package(Python 3 COMPONENTS Interpreter)
  set_package_properties(Python PROPERTIES TYPE OPTIONAL PURPOSE "Used to build documentation files")
  if(Python_EXECUTABLE AND VOTCA_XML_PARSER)
    file(GLOB XMLS ${CMAKE_CURRENT_SOURCE_DIR}/../xtp/xml/*.xml)

    set(XML_RST_FILES)
    foreach(_XML ${XMLS})
      get_filename_component(_FILE "${_XML}" NAME_WE)
      add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/xtp/${_FILE}.rst COMMAND ${Python_EXECUTABLE} ${VOTCA_XML_PARSER} -i ${_XML} -o ${VOTCA_SPHINX_DIR}/xtp/${_FILE}.rst
        DEPENDS ${VOTCA_XML_PARSER} ${_XML})
      list(APPEND XML_RST_FILES ${VOTCA_SPHINX_DIR}/xtp/${_FILE}.rst)
    endforeach()
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/xtp/qmpackage_defaults.rst COMMAND ${Python_EXECUTABLE} ${VOTCA_XML_PARSER} -m qm -i ${CMAKE_CURRENT_SOURCE_DIR}/../xtp/data/qmpackage_defaults.xml -o ${VOTCA_SPHINX_DIR}/xtp/qmpackage_defaults.rst
    DEPENDS ${VOTCA_XML_PARSER} ${CMAKE_CURRENT_SOURCE_DIR}/../xtp/data/qmpackage_defaults.xml)
    list(APPEND XML_RST_FILES ${VOTCA_SPHINX_DIR}/xtp/qmpackage_defaults.rst)
  endif()
  add_custom_target(doc-xtp-manual DEPENDS ${XML_RST_FILES} ${MANUAL_RST_FILES})
  add_dependencies(doc-xtp doc-xtp-manual)
endif()
