file(GLOB VOTCA_SCRIPTS *.sh *.pl *.py)

install(FILES ${VOTCA_FILES} DESTINATION ${DATA}/scripts/dipro)
install(PROGRAMS ${VOTCA_SCRIPTS} DESTINATION ${DATA}/scripts/dipro)
