foreach(PIC scheme schemes_all)
  #Run FIGTEX2EPS in ${PIC} subdir to avoid overwriting on latex tmp files in parallel build
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/${PIC})
  add_custom_command(OUTPUT ${CURRENT_DIR}/${PIC}.eps
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/${PIC}.fig
       ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/${PIC}/${PIC}.fig
    COMMAND ${FIGTEX2EPS} -v ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/${PIC}/${PIC}.fig
       ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/${PIC}/${PIC}.eps
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/${PIC}/${PIC}.eps
       ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/${PIC}.eps
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/${PIC}
    DEPENDS ${CURRENT_DIR}/${PIC}.fig)
endforeach()
