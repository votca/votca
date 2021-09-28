find_package(Python 3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

file(GLOB_RECURSE BIND_SOURCES binds/*.cc)

pybind11_add_module(xtp_binds ${BIND_SOURCES})
target_include_directories(xtp_binds PUBLIC binds/include)

target_link_libraries(xtp_binds PUBLIC VOTCA::votca_xtp)

if(NOT PYTHON_INSTDIR)
  execute_process(
    COMMAND
      ${Python_EXECUTABLE} -c
      "import distutils.sysconfig as cg; print(cg.get_python_lib(prefix='${CMAKE_INSTALL_PREFIX}', plat_specific=True, standard_lib=False))"
    OUTPUT_VARIABLE PYTHON_INSTDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(NOT PYTHON_INSTDIR)
install(TARGETS xtp_binds LIBRARY DESTINATION ${PYTHON_INSTDIR}/pyxtp)

file(GLOB_RECURSE PYXTP_FILES *.py)
install(FILES ${PYXTP_FILES} DESTINATION ${PYTHON_INSTDIR}/pyxtp)

if(ENABLE_TESTING)
  add_test(NAME integration_test_xtp_binds COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/examples.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_property(TEST integration_test_xtp_binds PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/..:$ENV{PYTHONPATH};VOTCASHARE=${CMAKE_CURRENT_SOURCE_DIR}/../../share")
  set_tests_properties(integration_test_xtp_binds PROPERTIES LABELS "csg;tools;votca;xtp;integration")
endif()
