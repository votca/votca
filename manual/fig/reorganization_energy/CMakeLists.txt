foreach(PIC marcus_parabolas monomer_parabolas) 
  add_custom_command(OUTPUT ${CURRENT_DIR}/${PIC}.eps
    COMMAND ${INKSCAPE} -f ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/${PIC}.svg
      -E ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/${PIC}.eps
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}
    DEPENDS ${CURRENT_DIR}/${PIC}.svg)
endforeach()
