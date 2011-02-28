file(GLOB CSG_BO_SOURCES *.cc)
add_executable(csg_boltzmann ${CSG_BO_SOURCES})
target_link_libraries(csg_boltzmann votca_csg)
install(TARGETS csg_boltzmann RUNTIME DESTINATION bin)
