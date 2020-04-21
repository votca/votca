find_package(Doxygen)
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
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
  file(GLOB PAGES *.md *.rst)
  file(GLOB IMAGES *.png)	
  configure_file(conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/sphinx/conf.py @ONLY)
  foreach(_MODULE ${ENABLED_VOTCA_PACKAGES})
    file(GLOB _MODULE_PAGES ${PROJECT_SOURCE_DIR}/${_MODULE}/guide/*.md ${PROJECT_SOURCE_DIR}/${_MODULE}/share/doc/*.md)
    list(APPEND PAGES ${_MODULE_PAGES})
  endforeach()
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/sphinx/.pages
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PAGES} ${CMAKE_CURRENT_BINARY_DIR}/sphinx
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${IMAGES} ${CMAKE_CURRENT_BINARY_DIR}/sphinx



    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/sphinx/.pages
    DEPENDS ${PAGES}
  )
  add_custom_target(doc COMMAND ${SPHINX_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/sphinx ${CMAKE_CURRENT_BINARY_DIR}/sphinx.html
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/sphinx/.pages ${CMAKE_CURRENT_BINARY_DIR}/sphinx/conf.py)
endif()
