foreach(PIC flow_boltzmann flow_fmatch flow_ibi)
  add_custom_command(OUTPUT ${CURRENT_DIR}/${PIC}.eps
    COMMAND ${TEX2EPS} ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/${PIC}.tex
      ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/${PIC}.eps
    COMMAND ${BBOX_ADD} ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/${PIC}.eps
    DEPENDS ${CURRENT_DIR}/${PIC}.tex)
endforeach()
