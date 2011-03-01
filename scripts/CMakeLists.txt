option(WITH_RC_FILES "Install votca rc files, no need when installing under /usr" ON)
if (WITH_RC_FILES)
  configure_file(VOTCARC.csh.in ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.csh @ONLY)
  configure_file(VOTCARC.bash.in ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.bash @ONLY)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.csh
    ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.bash DESTINATION bin)
endif (WITH_RC_FILES)
