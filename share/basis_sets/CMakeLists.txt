file(GLOB_RECURSE VOTCA_XTP_XML *.xml)
install(FILES ${VOTCA_XTP_XML} DESTINATION ${DATA}/xtp/basis_sets)
