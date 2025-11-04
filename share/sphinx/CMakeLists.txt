if(NOT SPHINX_FOUND)
  return()
endif()

set(SPHINX_DEPS ${VOTCA_SPHINX_DIR}/conf.py)
file(GLOB PAGES *.rst)
file(GLOB IMAGES *.png)
configure_file(conf.py.in ${VOTCA_SPHINX_DIR}/conf.py @ONLY)
foreach(_DEP ${PAGES} ${IMAGES})
  get_filename_component(_FILE "${_DEP}" NAME)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/${_FILE}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_DEP} ${VOTCA_SPHINX_DIR}/${_FILE} DEPENDS ${_DEP}
  )
  list(APPEND SPHINX_DEPS ${VOTCA_SPHINX_DIR}/${_FILE})
endforeach()
file(COPY _themes DESTINATION ${VOTCA_SPHINX_DIR})
foreach(_MODULE csg xtp xtp-tutorials)
  file(MAKE_DIRECTORY ${VOTCA_SPHINX_DIR}/${_MODULE})
  file(GLOB ${_MODULE}_FILES ${_MODULE}/*)
  add_custom_target(${_MODULE}_copy_rst
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/${_MODULE} ${VOTCA_SPHINX_DIR}/${_MODULE}
    DEPENDS ${${_MODULE}_FILES})
    add_dependencies(sphinx-deps ${_MODULE}_copy_rst)
endforeach()
add_custom_target(sphinx COMMAND ${SPHINX_EXECUTABLE} ${VOTCA_SPHINX_DIR} ${VOTCA_SPHINX_OUTPUT_DIR}
  DEPENDS ${SPHINX_DEPS})
add_dependencies(sphinx sphinx-deps)
