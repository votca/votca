find_package(Doxygen)
set_package_properties(Doxygen PROPERTIES TYPE OPTIONAL PURPOSE "Needed to build API documentation")
if (DOXYGEN_FOUND)
  set(DOXYGEN_EXTRACT_ALL YES)
  set(DOXYGEN_EXTRACT_PRIVATE YES)
  set(DOXYGEN_EXTRACT_PRIV_VIRTUAL YES)
  set(DOXYGEN_EXTRACT_STATIC YES)
  set(DOXYGEN_EXCLUDE_PATTERNS "*.py")
  set(DOXYGEN_SOURCE_BROWSER YES)
  # when using cmake-3.16 add USE_STAMP_FILE below
  doxygen_add_docs(doxygen ${ENABLED_VOTCA_PACKAGES}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Build doxygen documentation")
endif()

find_program(SPHINX_EXECUTABLE NAMES sphinx-build DOC "Sphinx documentation generation tool (http://www.sphinx-doc.org/)")
find_package_handle_standard_args(SPHINX REQUIRED_VARS SPHINX_EXECUTABLE)
if(SPHINX_FOUND)
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
  file(COPY _static DESTINATION ${VOTCA_SPHINX_DIR})
  add_custom_target(doc COMMAND ${SPHINX_EXECUTABLE} ${VOTCA_SPHINX_DIR} ${VOTCA_SPHINX_OUTPUT_DIR}
    DEPENDS ${SPHINX_DEPS})
  foreach(_MODULE ${ENABLED_VOTCA_PACKAGES})
    if(TARGET doc-${_MODULE})
      add_dependencies(doc doc-${_MODULE})
    endif()
  endforeach()
endif()
