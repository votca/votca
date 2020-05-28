
if(NOT VOTCA_SPHINX_DIR)
  set(VOTCA_SPHINX_DIR "${CMAKE_CURRENT_BINARY_DIR}")
endif()

file(GLOB NOTEBOOKS "${PROJECT_SOURCE_DIR}/GROMACS/Methane/QMMM_GROMACS.ipynb")

set(IPYNB_FILES)
foreach(_DEP ${NOTEBOOKS})
  get_filename_component(_FILE "${_DEP}" NAME)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/${_FILE}
    COMMAND ${CMAKE_COMMAND} -E copy ${_DEP} ${VOTCA_SPHINX_DIR}/${_FILE}
    )
  list(APPEND IPYNB_FILES ${VOTCA_SPHINX_DIR}/${_FILE})
endforeach()

message(STATUS "IPYNB_FILES: ${IPYNB_FILES}")
add_custom_target(doc-xtp-tutorials DEPENDS ${IPYNB_FILES})
