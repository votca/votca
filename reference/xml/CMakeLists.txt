add_custom_command(OUTPUT ${CURRENT_DIR}/config.t2t
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/../config.t2t ${CURRENT_DIR}/config.t2t
  DEPENDS ${CURRENT_DIR}/../config.t2t)

foreach(XML mapping.xml topol.xml csg_defaults.xml)
  add_custom_command(OUTPUT ${CURRENT_DIR}/${XML}.tex
    COMMAND CSG_PROPERTY=${CSG_PROPERTY} ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/xml2t2t.sh 
      ${CSG_SHARE}/xml/${XML} > ${CURRENT_DIR}/${XML}.t2t
    COMMAND ${TXT2TAGS_EXECUTABLE} -t tex --no-header --infile ${CURRENT_DIR}/${XML}.t2t
      --outfile ${CURRENT_DIR}/${XML}.tex
    DEPENDS ${CSG_SHARE}/xml/${XML} ${CURRENT_DIR}/config.t2t)
endforeach()
