foreach(PROG test_elements)
  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)
  add_executable(unit_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(unit_${PROG} votca_tools)

  if(ENABLE_TESTING)
    add_test(run_unit_${PROG} unit_${PROG})
    # run tests for tools (for coverage) as well
    set_tests_properties(run_unit_${PROG} PROPERTIES LABELS "tools;votca")
  endif(ENABLE_TESTING)
endforeach(PROG)
