add_custom_command(OUTPUT ${CURRENT_DIR}/config.t2t
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/../config.t2t ${CURRENT_DIR}/config.t2t
  DEPENDS ${CURRENT_DIR}/../config.t2t)

set(VOTCASHARE ${CMAKE_CURRENT_SOURCE_DIR}/../share)
foreach(XML map.xml options.xml)
  add_custom_command(OUTPUT ${CURRENT_DIR}/${XML}.tex
    COMMAND CSG_PROPERTY=${CSG_PROPERTY} ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/xml2t2t.sh 
      ${VOTCASHARE}/xtp/xml/${XML} > ${CURRENT_DIR}/${XML}.t2t
    COMMAND ${TXT2TAGS_EXECUTABLE} -t tex --no-header --infile ${CURRENT_DIR}/${XML}.t2t
      --outfile ${CURRENT_DIR}/${XML}.tex
    DEPENDS ${VOTCASHARE}/xtp/xml/${XML} ${CURRENT_DIR}/config.t2t)
endforeach()
