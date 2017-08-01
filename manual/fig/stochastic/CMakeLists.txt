add_custom_command(OUTPUT ${CURRENT_DIR}/overview_stochastic.eps
  COMMAND ${INKSCAPE} -f ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/overview_stochastic.svg
    -E ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/overview_stochastic.eps
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}
  DEPENDS ${CURRENT_DIR}/overview_stochastic.svg)
