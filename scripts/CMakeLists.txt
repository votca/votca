set(XTP_RST_FILES)
set(XTP_SCRIPTS)
find_package(Python 3 COMPONENTS Interpreter)
if(NOT Python_EXECUTABLE)
  return()
endif()
foreach(_pkg lxml h5py rdkit)
  if(NOT PY${_pkg}_FOUND)
    execute_process(COMMAND ${Python_EXECUTABLE} -c "import ${_pkg}" RESULT_VARIABLE IMPORT_py${_pkg})
    if(IMPORT_py${_pkg} EQUAL 0)
      set(PY${_pkg}_FOUND SUCCESS TRUE CACHE INTERNAL "")
      message(STATUS "Found python module ${_pkg}")
    endif()
  endif()
endforeach()
if(PYxml_FOUND)
  list(APPEND XTP_SCRIPTS xtp_basisset xtp_update_mapfile)
endif()
if(PYrdkit_FOUND)
  list(APPEND XTP_SCRIPTS xtp_autogen_mapping)
endif()
if(PYh5py_FOUND)
  list(APPEND XTP_SCRIPTS xtp_qmmm2qm)
endif()
foreach(SCRIPT ${XTP_SCRIPTS}) 
  configure_file(${SCRIPT}.in ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.in2 @ONLY)
  add_custom_target(${SCRIPT}_build ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT})
  add_custom_command(OUTPUT ${SCRIPT} COMMAND ${CMAKE_COMMAND}
    -DINPUT="${SCRIPT}.in2" -DOUTPUT="${SCRIPT}"
    -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
    -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${PROJECT_SOURCE_DIR}/CMakeModules/gitscript.cmake
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.in2 ${PROJECT_SOURCE_DIR}/CMakeModules/gitscript.cmake)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${SCRIPT})
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} DESTINATION ${CMAKE_INSTALL_BINDIR})
  if (BUILD_MANPAGES)
    add_custom_command(OUTPUT ${SCRIPT}.man
      COMMAND $<TARGET_FILE:VOTCA::votca_help2doc> --name ${SCRIPT} --format groff --out ${SCRIPT}.man
      COMMENT "Building ${SCRIPT} manpage"
      DEPENDS $<TARGET_FILE:VOTCA::votca_help2doc> ${SCRIPT}_build)
    add_custom_target(${SCRIPT}_manpage DEPENDS ${SCRIPT}.man)
    add_dependencies(manpages ${SCRIPT}_manpage)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME ${SCRIPT}.1)
  endif()
  if(VOTCA_SPHINX_DIR)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/xtp/${SCRIPT}.rst
      COMMAND $<TARGET_FILE:VOTCA::votca_help2doc> --name ${SCRIPT} --format rst --out ${VOTCA_SPHINX_DIR}/xtp/${SCRIPT}.rst
      COMMENT "Building ${SCRIPT} rst doc"
      DEPENDS $<TARGET_FILE:VOTCA::votca_help2doc> ${SCRIPT})
    list(APPEND TOOLS_RST_FILES ${VOTCA_SPHINX_DIR}/xtp/${SCRIPT}.rst)
  endif(VOTCA_SPHINX_DIR)
endforeach(SCRIPT)

if(VOTCA_SPHINX_DIR)
  file(MAKE_DIRECTORY ${VOTCA_SPHINX_DIR}/xtp)
  add_custom_target(doc-xtp-scripts DEPENDS ${TOOLS_RST_FILES})
  add_dependencies(doc-xtp doc-xtp-scripts)
endif()
