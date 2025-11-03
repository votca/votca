get_property(VOTCA_BINARIES TARGET votca_csg PROPERTY BINARIES)
list(SORT VOTCA_BINARIES)
# for .bash we need a space separated list
list(JOIN VOTCA_BINARIES " " VOTCA_BINARIES)
configure_file(csg-completion.bash.in "${CMAKE_CURRENT_BINARY_DIR}/csg-completion.bash" @ONLY)
if(INSTALL_RC_FILES)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/csg-completion.bash DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/rc)
endif()

