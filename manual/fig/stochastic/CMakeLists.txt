add_custom_command(OUTPUT ${CURRENT_DIR}/overview_stochastic.eps
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/overview_stochastic.svg
    ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/overview_stochastic.svg
  COMMAND ${INKSCAPE} -f overview_stochastic.svg -E overview_stochastic.eps
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}
  DEPENDS ${CURRENT_DIR}/overview_stochastic.svg)
