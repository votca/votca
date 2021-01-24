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
      COMMAND $<TARGET_FILE:VOTCA::votca_help2doc> --name $<TARGET_FILE:VOTCA::${SCRIPT}> --format groff --out ${SCRIPT}.man
      COMMENT "Building ${SCRIPT} manpage"
      DEPENDS $<TARGET_FILE:VOTCA::votca_help2doc> ${SCRIPT})
    add_custom_target(${SCRIPT}_manpage DEPENDS ${SCRIPT}.man)
    add_dependencies(manpages ${SCRIPT}_manpage)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME ${SCRIPT}.1)
  endif()
  if(VOTCA_SPHINX_DIR)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/csg/${SCRIPT}.rst
      COMMAND $<TARGET_FILE:VOTCA::votca_help2doc> --name $<TARGET_FILE:VOTCA::${SCRIPT}> --format rst --out ${VOTCA_SPHINX_DIR}/csg/${SCRIPT}.rst
      COMMENT "Building ${SCRIPT} rst doc"
      DEPENDS $<TARGET_FILE:VOTCA::votca_help2doc> ${SCRIPT})
    list(APPEND CSG_RST_FILES ${VOTCA_SPHINX_DIR}/csg/${SCRIPT}.rst)
  endif(VOTCA_SPHINX_DIR)
  add_executable(VOTCA::${SCRIPT} IMPORTED GLOBAL)
  set_property(TARGET VOTCA::${SCRIPT} PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}")
endforeach(SCRIPT)

set_target_properties(VOTCA::csg_call PROPERTIES VOTCASHARE ${PROJECT_SOURCE_DIR}/share)
set_target_properties(VOTCA::csg_call PROPERTIES VOTCA_PATH "${CMAKE_CURRENT_BINARY_DIR}:${PROJECT_BINARY_DIR}/src/tools")
set_target_properties(VOTCA::csg_inverse PROPERTIES VOTCA_CSG_DEFAULTS "${PROJECT_BINARY_DIR}/share/xml/csg_defaults.xml")

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/csg-completion.bash DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/rc)

if(VOTCA_SPHINX_DIR)
  add_custom_target(doc-csg-rst-scripts DEPENDS ${CSG_RST_FILES})
  add_dependencies(doc-csg doc-csg-rst-scripts)
endif()
