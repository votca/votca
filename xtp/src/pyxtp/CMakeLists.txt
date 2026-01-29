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

#find the user site-packages directory
execute_process(
  COMMAND ${Python_EXECUTABLE} -c "import site; print(site.getusersitepackages())"
  OUTPUT_VARIABLE PYTHON_USER_SITE
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Python user site: ${PYTHON_USER_SITE}")

#installation directory is the user's site-packages directory
if(NOT CMAKE_INSTALL_PYTHON_LIBDIR)
  set(CMAKE_INSTALL_PYTHON_LIBDIR "${PYTHON_USER_SITE}")
endif()


install(TARGETS xtp_binds LIBRARY DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}/pyxtp_binds)
add_library(VOTCA::xtp_binds ALIAS xtp_binds)

install (CODE "
    message(STATUS \"Installing Python package with pip from ${CMAKE_CURRENT_SOURCE_DIR} to ${CMAKE_INSTALL_PYTHON_LIBDIR}...\")
    execute_process(
    COMMAND ${Python_EXECUTABLE} -m pip install --upgrade --target ${CMAKE_INSTALL_PYTHON_LIBDIR} ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE _pip_result
  )
  if(NOT _pip_result EQUAL 0)
    message(FATAL_ERROR \"pip install failed with code \${_pip_result}\")
  endif()"
)

# copy pyxtp to build dir as well.
add_custom_target(pyxtp_bindir_copy COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pyxtp ${CMAKE_CURRENT_BINARY_DIR})
add_dependencies(xtp_binds pyxtp_bindir_copy)
set_target_properties(xtp_binds PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${XTP_PYPATH}/pyxtp_binds)


if(BUILD_TESTING AND PYpytest_FOUND)
  add_test(NAME integration_test_xtp_binds COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR} )
  set_property(TEST integration_test_xtp_binds PROPERTY ENVIRONMENT "PYTHONPATH=${XTP_PYPATH}:$ENV{PYTHONPATH};VOTCASHARE=${CMAKE_CURRENT_SOURCE_DIR}/../../share")
endif()
