file(GLOB_RECURSE VOTCA_HEADERS *.h potentialfunctions/*.h)
install(FILES ${VOTCA_HEADERS} DESTINATION include/votca/csg)
