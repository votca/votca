set(CSG_DEFAULTS csg_defaults.xml)
set(MDRUN "gmx mdrun")
set(G_ENERGY "gmx energy")
set(GROMPP "gmx grompp")

configure_file(${CSG_DEFAULTS}.in ${CMAKE_CURRENT_BINARY_DIR}/${CSG_DEFAULTS} @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${CSG_DEFAULTS} DESTINATION ${DATA}/xml)

file(GLOB_RECURSE VOTCA_XML mapping.xml)
install(FILES ${VOTCA_XML} DESTINATION ${DATA}/xml)

file(GLOB_RECURSE TOPOL_XML topol.xml)
install(FILES ${TOPOL_XML} DESTINATION ${DATA}/xml)
