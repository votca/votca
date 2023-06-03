find_package(Python 3.5 COMPONENTS Development)
set_package_properties(Python PROPERTIES TYPE RECOMMENDED PURPOSE "Used to build xtp python bindings")
find_package(pybind11)
set_package_properties(pybind11 PROPERTIES TYPE RECOMMENDED PURPOSE "Used to build xtp python bindings")

foreach(_pkg ase xmltodict pytest)
  if(NOT PY${_pkg}_FOUND)
    execute_process(COMMAND ${Python_EXECUTABLE} -c "import ${_pkg}" RESULT_VARIABLE IMPORT_py${_pkg})
    if(IMPORT_py${_pkg} EQUAL 0)
      set(PY${_pkg}_FOUND TRUE CACHE INTERNAL "")
      message(STATUS "Found python module ${_pkg}")
    endif()
  endif()
endforeach()

if(NOT pybind11_FOUND OR NOT Python_Development_FOUND OR NOT PYase_FOUND OR NOT PYxmltodict_FOUND)
  message(STATUS "Skipping xtp python bindings")
  return()
endif()

file(GLOB_RECURSE BIND_SOURCES binds/*.cc)

pybind11_add_module(xtp_binds ${BIND_SOURCES})
target_include_directories(xtp_binds PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/binds/include>)

target_link_libraries(xtp_binds PUBLIC VOTCA::votca_xtp)

if(NOT CMAKE_INSTALL_PYTHON_LIBDIR)
  set(CMAKE_INSTALL_PYTHON_LIBDIR "${CMAKE_INSTALL_LIBDIR}/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages")
endif()
install(TARGETS xtp_binds LIBRARY DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}/pyxtp)

file(GLOB_RECURSE PYXTP_FILES pyxtp/*.py)
message(STATUS "PYTHON INSTALL: ${CMAKE_INSTALL_PYTHON_LIBDIR}")
install(FILES ${PYXTP_FILES} DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}/pyxtp)

# copy py files to build dir, so CMAKE_CURRENT_BINARY_DIR is a complete py module with script and xtp_binds
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

if(BUILD_TESTING AND PYpytest_FOUND)
  add_test(NAME integration_test_xtp_binds COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR} )
  set_property(TEST integration_test_xtp_binds PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/../:$ENV{PYTHONPATH};VOTCASHARE=${CMAKE_CURRENT_SOURCE_DIR}/../../share")
endif()
