add_subdirectory(libxtp)
add_subdirectory(tools)
add_subdirectory(tests)

if(TARGET Python::Module AND pybind11_FOUND)
  add_subdirectory(pyxtp)
endif()
