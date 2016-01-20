configure_file(votca_xtp_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_xtp_config.h)

file(GLOB_RECURSE VOTCA_HEADERS *.h ${CMAKE_CURRENT_BINARY_DIR}/votca_xtp_config.h)
install(FILES ${VOTCA_HEADERS} DESTINATION include/votca/xtp)
