option(INSTALL_RC_FILES "Install votca rc files, no need when installing under /usr" ON)
add_feature_info(INSTALL_RC_FILES INSTALL_RC_FILES "Install votca rc files (VOTCARC.bash etc.)")
if(INSTALL_RC_FILES)
  if(APPLE)
    set(LD_LIBRARY_PATH "DYLD_LIBRARY_PATH")
  else(APPLE)
    set(LD_LIBRARY_PATH "LD_LIBRARY_PATH")
  endif(APPLE)
  configure_file(VOTCARC.bash.in ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.bash @ONLY)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.bash DESTINATION ${CMAKE_INSTALL_BINDIR})
endif(INSTALL_RC_FILES)

set(TOOLS_RST_FILES)
foreach(SCRIPT votca_compare votca_help2doc)
  configure_file(${SCRIPT}.in ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.tmp.out @ONLY)
  add_custom_target(${SCRIPT}_build ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT})
  add_custom_command(OUTPUT ${SCRIPT} COMMAND ${CMAKE_COMMAND}
    -DINPUT="${SCRIPT}.tmp.out" -DOUTPUT="${SCRIPT}"
    -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
    -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${PROJECT_SOURCE_DIR}/CMakeModules/gitscript.cmake
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.tmp.out ${PROJECT_SOURCE_DIR}/CMakeModules/gitscript.cmake)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${SCRIPT})
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(VOTCA::${SCRIPT} IMPORTED GLOBAL)
  set_property(TARGET VOTCA::${SCRIPT} PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}")
  add_dependencies(VOTCA::${SCRIPT} ${SCRIPT}_build)
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
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/tools/${SCRIPT}.rst
      COMMAND $<TARGET_FILE:VOTCA::votca_help2doc> --name $<TARGET_FILE:VOTCA::${SCRIPT}> --format rst --out ${VOTCA_SPHINX_DIR}/tools/${SCRIPT}.rst
      COMMENT "Building ${SCRIPT} rst doc"
      DEPENDS $<TARGET_FILE:VOTCA::votca_help2doc> ${SCRIPT})
    list(APPEND TOOLS_RST_FILES ${VOTCA_SPHINX_DIR}/tools/${SCRIPT}.rst)
  endif(VOTCA_SPHINX_DIR)
endforeach(SCRIPT)

if(VOTCA_SPHINX_DIR)
  add_custom_target(doc-tools-scripts DEPENDS ${TOOLS_RST_FILES})
  add_dependencies(doc-tools doc-tools-scripts)
endif()
