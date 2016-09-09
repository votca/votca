configure_file(votca_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)

if (NOT GSL_FOUND)
    file(GLOB_RECURSE NOT_VOTCA_HEADERS votca_gsl_boost_ublas_matrix_prod.h)
    list(REMOVE_ITEM VOTCA_HEADERS ${NOT_VOTCA_HEADERS})
endif (NOT GSL_FOUND)

file(GLOB_RECURSE VOTCA_HEADERS *.h ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)
install(FILES ${VOTCA_HEADERS} DESTINATION include/votca/xtp)
