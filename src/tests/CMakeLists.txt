
find_package(Boost 1.39.0 REQUIRED COMPONENTS unit_test_framework)
foreach(PROG test_glink test_aomatrix test_orbitals test_convergenceacc test_adiis test_diis test_atom test_eigen test_numericalintegration test_symmetric_matrix test_aobasis test_threecenter)
	file(GLOB ${PROG}_SOURCES ${PROG}*.cc)
	add_executable(unit_${PROG} ${${PROG}_SOURCES})
	target_link_libraries(unit_${PROG} votca_xtp ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${MKL_LIBRARIES})
	target_compile_definitions(unit_${PROG} PRIVATE BOOST_TEST_DYN_LINK)
	add_test(unit_${PROG} unit_${PROG})
	# run tests for xtp (for coverage) as well
	set_tests_properties(unit_${PROG} PROPERTIES LABELS "xtp;votca")
endforeach(PROG)

