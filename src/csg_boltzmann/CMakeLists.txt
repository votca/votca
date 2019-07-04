file(GLOB CSG_BO_SOURCES *.cc)
add_executable(csg_boltzmann ${CSG_BO_SOURCES})
target_link_libraries(csg_boltzmann votca_csg)
install(TARGETS csg_boltzmann EXPORT VOTCA_CSG_Targets RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
add_executable(VOTCA::csg_boltzmann ALIAS csg_boltzmann)
if (TXT2TAGS_FOUND AND BASH)
  add_custom_command(OUTPUT csg_boltzmann.man 
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/csg_boltzmann --help > csg_boltzmann.help DEPENDS csg_boltzmann
    COMMAND ${BASH} ${CMAKE_CURRENT_BINARY_DIR}/../../scripts/help2t2t csg_boltzmann.help > csg_boltzmann.t2t
    COMMAND ${TXT2TAGS_EXECUTABLE} -q -t man -i csg_boltzmann.t2t -o csg_boltzmann.man
    DEPENDS help2t2t_build csg_boltzmann)
  add_custom_target(csg_boltzmann_manpage DEPENDS csg_boltzmann.man)
  add_dependencies(manpages csg_boltzmann_manpage)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/csg_boltzmann.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME csg_boltzmann.1)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES csg_boltzmann.help csg_boltzmann.t2t)
endif(TXT2TAGS_FOUND AND BASH)

list(REMOVE_ITEM CSG_BO_SOURCES "main.cc")
add_library(votca_csg_boltzmann ${CSG_BO_SOURCES})
target_link_libraries(votca_csg_boltzmann votca_csg)
