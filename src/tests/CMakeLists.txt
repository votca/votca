find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)

foreach(PROG
  test_bead
  test_beadtriple
  test_basebead
  test_beadmotif_algorithms
  test_beadmotif_base
  test_beadmotif_connector
  test_beadstructure_base
  test_beadstructure_algorithms
  test_bondedstatistics
  test_csg_topology
  test_interaction
  test_lammpsdatareader 
  test_lammpsdumpreaderwriter
  test_nblist_3body
  test_nblistgrid_3body
  test_boundarycondition
  test_pdbreader
  test_tabulatedpotential
  test_triplelist )

  file(GLOB ${PROG}_SOURCES ${PROG}.cc)
  add_executable(unit_${PROG} ${${PROG}_SOURCES})
  if(PROG STREQUAL "test_bondedstatistics" OR PROG STREQUAL "test_tabulatedpotential")
    target_link_libraries(unit_${PROG} votca_csg_boltzmann)
  endif()
  target_link_libraries(unit_${PROG} votca_csg Boost::unit_test_framework)
  target_compile_definitions(unit_${PROG} PRIVATE BOOST_TEST_DYN_LINK)
  add_test(unit_${PROG} unit_${PROG})
  # run tests for csg (for coverage) as well
  set_tests_properties(unit_${PROG} PROPERTIES LABELS "csg;votca;unit")
endforeach(PROG)
