option(INSTALL_RC_FILES "Install votca rc files, no need when installing under /usr" ON)
add_feature_info(INSTALL_RC_FILES INSTALL_RC_FILES "Install votca rc files (VOTCARC.bash etc.)")
if(INSTALL_RC_FILES)
  if(APPLE)
    set(LD_LIBRARY_PATH "DYLD_LIBRARY_PATH")
  else(APPLE)
    set(LD_LIBRARY_PATH "LD_LIBRARY_PATH")
  endif(APPLE)
  configure_file(VOTCARC.csh.in ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.csh @ONLY)
  configure_file(VOTCARC.bash.in ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.bash @ONLY)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.csh
    ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.bash DESTINATION ${CMAKE_INSTALL_BINDIR})
endif(INSTALL_RC_FILES)

foreach(SCRIPT votca_compare)
  configure_file(${SCRIPT}.in ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.tmp.out @ONLY)
  add_custom_target(${SCRIPT}_build ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT})
  add_custom_command(OUTPUT ${SCRIPT} COMMAND ${CMAKE_COMMAND}
    -DINPUT="${SCRIPT}.tmp.out" -DOUTPUT="${SCRIPT}"
    -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
    -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${PROJECT_SOURCE_DIR}/CMakeModules/gitscript.cmake
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.tmp.out ${PROJECT_SOURCE_DIR}/CMakeModules/gitscript.cmake)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${SCRIPT})
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach(SCRIPT)

add_executable(VOTCA::votca_compare IMPORTED GLOBAL)
set_property(TARGET VOTCA::votca_compare PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/votca_compare")
