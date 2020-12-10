if(NOT VOTCA_SPHINX_DIR)
  set(VOTCA_SPHINX_DIR "${CMAKE_CURRENT_BINARY_DIR}")
endif()

file(GLOB RSTS ${CMAKE_CURRENT_SOURCE_DIR}/*.rst)
file(GLOB PNGS ${CMAKE_CURRENT_SOURCE_DIR}/fig/*.png)
set(MANUAL_RST_FILES)
set(MANUAL_PNG_FILES)
foreach(_RST ${RSTS})
  get_filename_component(_FILE "${_RST}" NAME)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/csg/${_FILE}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_RST} ${VOTCA_SPHINX_DIR}/csg/${_FILE}
    DEPENDS ${_RST})
  list(APPEND MANUAL_RST_FILES ${VOTCA_SPHINX_DIR}/csg/${_FILE})
endforeach()
foreach(_PNG ${PNGS})
    get_filename_component(_FILE "${_PNG}" NAME)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/csg/fig/${_FILE}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_PNG} ${VOTCA_SPHINX_DIR}/csg/fig/${_FILE}
      DEPENDS ${_PNG})
  list(APPEND MANUAL_PNG_FILES ${VOTCA_SPHINX_DIR}/csg/fig/${_FILE})
endforeach()
add_custom_target(doc-csg-manul DEPENDS ${MANUAL_RST_FILES} ${MANUAL_PNG_FILES})
if(TARGET doc-csg)
  add_dependencies(doc-csg doc-csg-manul)
endif()
