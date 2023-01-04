file(GLOB VOTCA_XTP_XML *.xml)
install(FILES ${VOTCA_XTP_XML} DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/xtp/xml)
file(GLOB VOTCA_XTP_SUBPKG_XML subpackages/*.xml)
install(FILES ${VOTCA_XTP_SUBPKG_XML} DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/xtp/xml/subpackages)

if(NOT SPHINX_FOUND)
  return()
endif()
  
set(XML_RST_FILES)
foreach(_XML ${VOTCA_XTP_XML})
  get_filename_component(_FILE "${_XML}" NAME_WE)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/xtp/${_FILE}.rst
    COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::extract_xml_metadata> -i ${_XML} -o ${VOTCA_SPHINX_DIR}/xtp/${_FILE}.rst
    DEPENDS VOTCA::extract_xml_metadata ${_XML} ${VOTCA_XTP_SUBPKG_XML})
  list(APPEND XML_RST_FILES ${VOTCA_SPHINX_DIR}/xtp/${_FILE}.rst)
endforeach()
add_custom_target(xtp-xml-rst DEPENDS ${XML_RST_FILES})
add_dependencies(sphinx-deps xtp-xml-rst)
