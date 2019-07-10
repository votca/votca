foreach(PIC marcus_parabolas monomer_parabolas) 
  add_custom_command(OUTPUT ${CURRENT_DIR}/${PIC}.eps
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/${PIC}.svg
      ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/${PIC}.svg
    COMMAND ${INKSCAPE} -f ${PIC}.svg -E ${PIC}.eps
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}
    DEPENDS ${CURRENT_DIR}/${PIC}.svg)
endforeach()
