find_package(Python 3.5 COMPONENTS Development)
set_package_properties(Python PROPERTIES TYPE RECOMMENDED PURPOSE "Used to build xtp python bindings")
find_package(pybind11)
set_package_properties(pybind11 PROPERTIES TYPE RECOMMENDED PURPOSE "Used to build xtp python bindings")

foreach(_pkg ase xmltodict)
  if(NOT PY${_pkg}_FOUND)
    execute_process(COMMAND ${Python_EXECUTABLE} -c "import ${_pkg}" RESULT_VARIABLE IMPORT_py${_pkg})
    if(IMPORT_py${_pkg} EQUAL 0)
      set(PY${_pkg}_FOUND SUCCESS TRUE CACHE INTERNAL "")
      message(STATUS "Found python module ${_pkg}")
    endif()
  endif()
endforeach()

if(NOT pybind11_FOUND OR NOT Python_Development_FOUND OR NOT PYase_FOUND OR NOT PYxml2dict_FOUND)
  return()
endif()

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

foreach(_PYFILE ${PYXTP_FILES})
  get_filename_component(_FILE "${_PYFILE}" NAME)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_FILE}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_PYFILE} ${CMAKE_CURRENT_BINARY_DIR}/${_FILE}
    DEPENDS ${_DEP}
  )
  list(APPEND BINARY_SCRIPTS ${CMAKE_CURRENT_BINARY_DIR}/${_FILE})
endforeach()
add_custom_target(xtp_binary_scripts DEPENDS ${BINARY_SCRIPTS})
add_dependencies(xtp_binds xtp_binary_scripts)

if(BUILD_TESTING)
  add_test(NAME integration_test_xtp_binds COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR} )
  set_property(TEST integration_test_xtp_binds PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/../:$ENV{PYTHONPATH};VOTCASHARE=${CMAKE_CURRENT_SOURCE_DIR}/../../share")
endif()
