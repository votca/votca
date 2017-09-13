file(GLOB VOTCA_SCRIPTS *.sh *.pl *.py)
file(GLOB VOTCA_FILES csg_table *.pm *.m *.octave) 

install(FILES ${VOTCA_FILES} DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/scripts/inverse)
install(PROGRAMS ${VOTCA_SCRIPTS} DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/scripts/inverse)
