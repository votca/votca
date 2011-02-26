configure_file(VOTCARC.csh.in VOTCARC.csh @ONLY)
configure_file(VOTCARC.bash.in VOTCARC.bash @ONLY)
install(PROGRAMS VOTCARC.csh VOTCARC.bash DESTINATION bin)
