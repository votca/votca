set(VOTCASHARE "${PROJECT_SOURCE_DIR}/share")
set(VOTCA_PATH "${CMAKE_CURRENT_BINARY_DIR}:${PROJECT_BINARY_DIR}/src/tools")
set(VOTCA_CSG_DEFAULTS "${PROJECT_BINARY_DIR}/share/xml/csg_defaults.xml")

set(CSG_RST_FILES)
foreach(SCRIPT csg_call csg_inverse)
  configure_file(${SCRIPT}.in ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.in2 @ONLY)
  add_custom_target(${SCRIPT}_build ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT})
  add_custom_command(OUTPUT ${SCRIPT} COMMAND ${CMAKE_COMMAND}
    -DINPUT="${SCRIPT}.in2" -DOUTPUT="${SCRIPT}"
    -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
    -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${PROJECT_SOURCE_DIR}/CMakeModules/gitscript.cmake
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.in2 ${PROJECT_SOURCE_DIR}/CMakeModules/gitscript.cmake)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${SCRIPT} ${SCRIPT}.tmp ${SCRIPT}.tmp.x)
  set_property(TARGET votca_csg APPEND PROPERTY BINARIES "${SCRIPT}")
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} DESTINATION ${CMAKE_INSTALL_BINDIR})
  if (BUILD_MANPAGES)
    add_custom_command(OUTPUT ${SCRIPT}.man
      COMMAND ${CMAKE_COMMAND} -E env VOTCASHARE="${VOTCASHARE}" PATH="${VOTCA_PATH}:$ENV{PATH}" VOTCA_CSG_DEFAULTS="${VOTCA_CSG_DEFAULTS}"
        $<TARGET_FILE:VOTCA::votca_help2doc> --name $<TARGET_FILE:VOTCA::${SCRIPT}> --format groff --out ${SCRIPT}.man
      COMMENT "Building ${SCRIPT} manpage"
      DEPENDS $<TARGET_FILE:VOTCA::votca_help2doc> ${SCRIPT})
    add_custom_target(${SCRIPT}_manpage DEPENDS ${SCRIPT}.man)
    add_dependencies(manpages ${SCRIPT}_manpage)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME ${SCRIPT}.1)
  endif()
  if(VOTCA_SPHINX_DIR)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/csg/${SCRIPT}.rst
      COMMAND ${CMAKE_COMMAND} -E env VOTCASHARE="${VOTCASHARE}" PATH="${VOTCA_PATH}:$ENV{PATH}" VOTCA_CSG_DEFAULTS="${VOTCA_CSG_DEFAULTS}"
	$<TARGET_FILE:VOTCA::votca_help2doc> --name $<TARGET_FILE:VOTCA::${SCRIPT}> --format rst --out ${VOTCA_SPHINX_DIR}/csg/${SCRIPT}.rst
      COMMENT "Building ${SCRIPT} rst doc"
      DEPENDS $<TARGET_FILE:VOTCA::votca_help2doc> ${SCRIPT})
    list(APPEND CSG_RST_FILES ${VOTCA_SPHINX_DIR}/csg/${SCRIPT}.rst)
  endif(VOTCA_SPHINX_DIR)
  add_executable(VOTCA::${SCRIPT} IMPORTED GLOBAL)
  set_property(TARGET VOTCA::${SCRIPT} PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}")
endforeach(SCRIPT)

set_target_properties(VOTCA::csg_call PROPERTIES VOTCASHARE "${VOTCASHARE}")
set_target_properties(VOTCA::csg_call PROPERTIES VOTCA_PATH "${VOTCA_PATH}")
set_target_properties(VOTCA::csg_inverse PROPERTIES VOTCA_CSG_DEFAULTS "${VOTCA_CSG_DEFAULTS}")

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/csg-completion.bash DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/rc)

if(VOTCA_SPHINX_DIR)
  #don't call "csg_call -l" as csg_call isn't ready at cmake time
  file(READ ${VOTCASHARE}/scripts/inverse/csg_table CSG_SCRIPTS)
  string(REGEX REPLACE "#[^\n]*\n" "" CSG_SCRIPTS "${CSG_SCRIPTS}") # remove comments
  string(REGEX REPLACE " +" " " CSG_SCRIPTS "${CSG_SCRIPTS}") # multi-space in single space
  string(REGEX REPLACE "\n" ";" CSG_SCRIPTS "${CSG_SCRIPTS}") # turn into list
  string(REGEX REPLACE "; *;" ";" CSG_SCRIPTS "${CSG_SCRIPTS}") # remove empty line

  set(CSG_TABLE "${VOTCA_SPHINX_DIR}/csg/csg_table.rst")
  file(WRITE ${CSG_TABLE} ".. list-table::\n   :header-rows: 1\n   :align: center\n\n")
  file(APPEND ${CSG_TABLE} "   * - Key1\n     - Key2\n     - Scriptname\n")
  file(MAKE_DIRECTORY "${VOTCA_SPHINX_DIR}/csg/scripts")
  set(CSG_SCRIPTS_LIST)
  foreach(_LINE ${CSG_SCRIPTS})
    separate_arguments(_LINE)
    list(GET _LINE 0 TAG1)
    file(APPEND ${CSG_TABLE} "   * - ${TAG1}\n")
    list(GET _LINE 1 TAG2)
    file(APPEND ${CSG_TABLE} "     - ${TAG2}\n")
    list(GET _LINE 2 SCRIPT)
    list(REMOVE_AT _LINE 0 1)
    list(JOIN _LINE " " _LINE)
    file(APPEND ${CSG_TABLE} "     - ${_LINE}\n")
    set(FULL_SCRIPT "${VOTCA_SPHINX_DIR}/csg/scripts/${SCRIPT}.rst")
    if("${SCRIPT}" STREQUAL "")
      message("CJ: ${TAG1} ${TAG2} xx ${SCRIPT}")
    endif()
    list(FIND CSG_RST_FILES ${FULL_SCRIPT} _alreadyThere)
    if(NOT _alreadyThere STREQUAL "-1")
      continue()
    endif()
    list(APPEND CSG_RST_FILES ${FULL_SCRIPT})
    string(REGEX REPLACE "." "^" UNDERLINE "${SCRIPT}")
    list(APPEND CSG_SCRIPTS_LIST "${SCRIPT}\n${UNDERLINE}\n\n.. include:: scripts/${SCRIPT}.rst\n\n")
    add_custom_command(OUTPUT ${FULL_SCRIPT}
      COMMAND ${CMAKE_COMMAND} -E env VOTCASHARE="${VOTCASHARE}" PATH="${VOTCA_PATH}:$ENV{PATH}" VOTCA_CSG_DEFAULTS="${VOTCA_CSG_DEFAULTS}"
      $<TARGET_FILE:VOTCA::votca_help2doc> --name "$<TARGET_FILE:VOTCA::csg_call> ${TAG1} ${TAG2}" --format rst --out ${FULL_SCRIPT}
      COMMENT "Building ${SCRIPT} rst doc (from ${TAG1} ${TAG2})"
      DEPENDS $<TARGET_FILE:VOTCA::votca_help2doc> ${VOTCASHARE}/scripts/inverse/${SCRIPT})
  endforeach(_LINE)
  list(SORT CSG_SCRIPTS_LIST)
  set(CSG_SCRIPTS_RST "${VOTCA_SPHINX_DIR}/csg/csg_scripts.rst")
  file(WRITE ${CSG_SCRIPTS_RST} "\n")
  foreach(_LINE ${CSG_SCRIPTS_LIST})
    file(APPEND ${CSG_SCRIPTS_RST} ${_LINE})
  endforeach()

  add_custom_target(doc-csg-rst-scripts DEPENDS ${CSG_RST_FILES})
  add_dependencies(doc-csg doc-csg-rst-scripts)
endif()
