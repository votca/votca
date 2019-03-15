find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)

# Each test listed in Alphabetical order
foreach(PROG
		test_constants
		test_cubicspline
		test_datacollection
		test_edge_base
		test_edgecontainer
		test_elements
		test_floatingpointcomparison
		test_graph_base
		test_graph_bf_visitor
		test_graph_df_visitor
		test_graphalgorithm
		test_graphdistvisitor
		test_graphnode
		test_graphvisitor
		test_histogramnew
		test_identity
		test_linalg
		test_name
		test_reducededge
		test_reducedgraph
		test_table
		test_thread
		test_tokenizer
		test_typeconverter
        test_filesystem
		test_property)

  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)
  add_executable(unit_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(unit_${PROG} votca_tools ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  target_compile_definitions(unit_${PROG} PRIVATE BOOST_TEST_DYN_LINK)
  add_test(unit_${PROG} unit_${PROG})
  # run tests for tools (for coverage) as well
  set_tests_properties(unit_${PROG} PROPERTIES LABELS "tools;votca")

endforeach(PROG)
