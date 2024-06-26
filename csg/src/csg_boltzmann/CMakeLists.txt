file(GLOB CSG_BO_LIB_SOURCES *.cc)
file(GLOB CSG_BO_MAIN_SOURCES main.cc)
list(REMOVE_ITEM CSG_BO_LIB_SOURCES ${CSG_BO_MAIN_SOURCES})

add_library(votca_csg_boltzmann OBJECT ${CSG_BO_LIB_SOURCES})
target_link_libraries(votca_csg_boltzmann votca_csg)

add_executable(csg_boltzmann ${CSG_BO_MAIN_SOURCES})
target_link_libraries(csg_boltzmann votca_csg_boltzmann votca_csg)
install(TARGETS csg_boltzmann EXPORT VOTCA_CSG_Targets RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
set_property(TARGET votca_csg APPEND PROPERTY BINARIES "csg_boltzmann")
add_executable(VOTCA::csg_boltzmann ALIAS csg_boltzmann)

if (BUILD_MANPAGES)
  add_custom_command(OUTPUT csg_boltzmann.man
    COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_help2doc> --name $<TARGET_FILE:VOTCA::csg_boltzmann> --format groff --out csg_boltzmann.man
    COMMENT "Building csg_boltzmann manpage"
    DEPENDS VOTCA::votca_help2doc csg_boltzmann)
  add_custom_target(csg_boltzmann_manpage DEPENDS csg_boltzmann.man)
  add_dependencies(manpages csg_boltzmann_manpage)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/csg_boltzmann.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME csg_boltzmann.1)
endif()

if(SPHINX_FOUND)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/csg/csg_boltzmann.rst
    COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_help2doc> --name $<TARGET_FILE:VOTCA::csg_boltzmann> --format rst --out ${VOTCA_SPHINX_DIR}/csg/csg_boltzmann.rst
    COMMENT "Building csg_boltzmann rst doc"
    DEPENDS VOTCA::votca_help2doc csg_boltzmann)
  add_custom_target(csg-boltzmann-rst DEPENDS ${VOTCA_SPHINX_DIR}/csg/csg_boltzmann.rst)
  add_dependencies(sphinx-deps csg-boltzmann-rst)
endif()
