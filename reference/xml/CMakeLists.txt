add_custom_command(OUTPUT ${CURRENT_DIR}/config.t2t
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/../config.t2t ${CURRENT_DIR}/config.t2t
  DEPENDS ${CURRENT_DIR}/../config.t2t)

foreach(XML mapping.xml topol.xml)
  add_custom_command(OUTPUT ${CURRENT_DIR}/${XML}.tex
    COMMAND CSG_PROPERTY=$<TARGET_FILE:VOTCA::csg_property> ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/xml2t2t.sh
      ${CSG_SHARE}/xml/${XML} > ${CURRENT_DIR}/${XML}.t2t
    COMMAND ${TXT2TAGS_EXECUTABLE} -t tex --no-header --infile ${CURRENT_DIR}/${XML}.t2t
      --outfile ${CURRENT_DIR}/${XML}.tex
    DEPENDS ${CSG_SHARE}/xml/${XML} ${CURRENT_DIR}/config.t2t)
endforeach()
  
add_custom_command(OUTPUT ${CURRENT_DIR}/csg_defaults.xml.tex
  COMMAND CSG_PROPERTY=$<TARGET_FILE:VOTCA::csg_property> ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/xml2t2t.sh
    ${VOTCA_CSG_DEFAULTS} > ${CURRENT_DIR}/csg_defaults.xml.t2t
  COMMAND ${TXT2TAGS_EXECUTABLE} -t tex --no-header --infile
    ${CURRENT_DIR}/csg_defaults.xml.t2t --outfile ${CURRENT_DIR}/csg_defaults.xml.tex
  DEPENDS ${VOTCA_CSG_DEFAULTS} ${CURRENT_DIR}/config.t2t)
