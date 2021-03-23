file(GLOB VOTCA_SCRIPTS *.sh *.pl *.py)
file(GLOB VOTCA_FILES csg_table *.pm *.m *.octave)

install(FILES ${VOTCA_FILES} DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/scripts/inverse)
install(PROGRAMS ${VOTCA_SCRIPTS} DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/scripts/inverse)

if(ENABLE_TESTING)
  if(ENABLE_COVERAGE_BUILD)
    set(PY_COV_OPTS "--cov-report=xml --cov=.")
  endif()

  add_test(unit_iie.py ${Python_EXECUTABLE} -m pytest -q ${PY_COV_OPTS} ${CMAKE_CURRENT_SOURCE_DIR}/iie.py)
  set_tests_properties(unit_iie.py PROPERTIES LABELS "csg;votca;unit")
endif(ENABLE_TESTING)
