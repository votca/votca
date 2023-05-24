file(GLOB_RECURSE BIND_SOURCES binds/*.cc)

pybind11_add_module(xtp_binds ${BIND_SOURCES})
target_include_directories(xtp_binds PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/binds/include>)

target_link_libraries(xtp_binds PUBLIC VOTCA::votca_xtp)

if(NOT PYTHON_INSTDIR)
  execute_process(
    COMMAND
      ${Python_EXECUTABLE} -c
      "import distutils.sysconfig as cg; print(cg.get_python_lib(plat_specific=True, standard_lib=False))"
    OUTPUT_VARIABLE PYTHON_INSTDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(NOT PYTHON_INSTDIR)
install(TARGETS xtp_binds LIBRARY DESTINATION ${PYTHON_INSTDIR}/pyxtp)

file(GLOB_RECURSE PYXTP_FILES pyxtp/*.py)
message(STATUS "PYTHON INSTALL: ${PYTHON_INSTDIR}")
install(FILES ${PYXTP_FILES} DESTINATION ${PYTHON_INSTDIR}/pyxtp)

if(BUILD_TESTING)
  add_test(NAME integration_test_xtp_binds COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test/ WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/)
  set_property(TEST integration_test_xtp_binds PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}:$ENV{PYTHONPATH};VOTCASHARE=${CMAKE_CURRENT_SOURCE_DIR}/../../share")
  # list(APPEND PYTHONPATH "${PYTHON_INSTDIR}")
  # set_property(TEST integration_test_xtp_binds PROPERTY ENVIRONMENT "VOTCASHARE=${CMAKE_CURRENT_SOURCE_DIR}/../../share")
  set_tests_properties(integration_test_xtp_binds PROPERTIES LABELS "csg;tools;votca;xtp;integration")
endif()