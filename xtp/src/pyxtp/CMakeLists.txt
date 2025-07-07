find_package(Python 3.5 COMPONENTS Development Interpreter)
set_package_properties(Python PROPERTIES TYPE RECOMMENDED PURPOSE "Used to build xtp python bindings")
find_package(pybind11)
set_package_properties(pybind11 PROPERTIES TYPE RECOMMENDED PURPOSE "Used to build xtp python bindings")

if(NOT pybind11_FOUND OR NOT Python_Development_FOUND OR NOT Python_Interpreter_FOUND)
  message(STATUS "Skipping xtp python bindings")
  return()
endif()

file(GLOB_RECURSE BIND_SOURCES binds/*.cc)

pybind11_add_module(xtp_binds ${BIND_SOURCES})
target_include_directories(xtp_binds PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/binds/include>)
set(XTP_PYPATH "${CMAKE_CURRENT_BINARY_DIR}/../")
set_target_properties(xtp_binds PROPERTIES XTP_PYPATH "${XTP_PYPATH}")

target_link_libraries(xtp_binds PUBLIC VOTCA::votca_xtp)

if(NOT CMAKE_INSTALL_PYTHON_LIBDIR)
  set(CMAKE_INSTALL_PYTHON_LIBDIR "${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages")
endif()

install(TARGETS xtp_binds LIBRARY DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}/pyxtp_binds)

install (CODE "
    message(STATUS \"Installing Python package with pip from ${CMAKE_CURRENT_SOURCE_DIR}...\")
    execute_process(
    COMMAND ${Python_EXECUTABLE} -m pip install --upgrade --target ${CMAKE_INSTALL_PYTHON_LIBDIR} ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE _pip_result
  )
  if(NOT _pip_result EQUAL 0)
    message(FATAL_ERROR \"pip install failed with code \${_pip_result}\")
  endif()"
)

add_library(VOTCA::xtp_binds ALIAS xtp_binds)

if(BUILD_TESTING AND PYpytest_FOUND)
  add_test(NAME integration_test_xtp_binds COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR} )
  set_property(TEST integration_test_xtp_binds PROPERTY ENVIRONMENT "PYTHONPATH=${XTP_PYPATH}:$ENV{PYTHONPATH};VOTCASHARE=${CMAKE_CURRENT_SOURCE_DIR}/../../share")
endif()
