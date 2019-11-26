find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)

# Each test listed in Alphabetical order
foreach(PROG
test_kspace
)

  file(GLOB ${PROG}_SOURCES ${PROG}*.cpp)
  add_executable(unit_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(unit_${PROG} ${QDMEwald_LIBS} Boost::unit_test_framework)
  target_compile_definitions(unit_${PROG} PRIVATE BOOST_TEST_DYN_LINK)
  add_test(unit_${PROG} unit_${PROG})
  set_tests_properties(unit_${PROG} PROPERTIES LABELS "qdmewald;unit")
endforeach(PROG)
