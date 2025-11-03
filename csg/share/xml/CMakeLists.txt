file(GLOB_RECURSE VOTCA_XML *.xml)
install(FILES ${VOTCA_XML} DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/xml)

if(SPHINX_FOUND)
  set(XML_RST_FILES)
  foreach(_XML ${VOTCA_XML})
    get_filename_component(_FILE "${_XML}" NAME_WE)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/csg/${_FILE}.rst COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::extract_xml_metadata> -i ${_XML} -m csg -o ${VOTCA_SPHINX_DIR}/csg/${_FILE}.rst
      DEPENDS VOTCA::extract_xml_metadata ${_XML})
    list(APPEND XML_RST_FILES ${VOTCA_SPHINX_DIR}/csg/${_FILE}.rst)
  endforeach()
  add_custom_target(csg-xml-rst DEPENDS ${XML_RST_FILES})
  add_dependencies(sphinx-deps csg-xml-rst)
endif()
