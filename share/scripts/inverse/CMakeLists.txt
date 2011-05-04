file(GLOB VOTCA_SCRIPTS *.sh *.pl *.py)
file(GLOB VOTCA_FILES csg_table *.pm *.m *.octave) 

install(FILES ${VOTCA_FILES} DESTINATION @DATA@/scripts/inverse)
install(PROGRAMS ${VOTCA_SCRIPTS} DESTINATION @DATA@/scripts/inverse)
