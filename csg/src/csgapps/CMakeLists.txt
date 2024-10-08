option(INSTALL_CSGAPPS "Build csgapps" OFF)
add_feature_info(INSTALL_CSGAPPS INSTALL_CSGAPPS "Install extra app using csg")

file(GLOB FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} */*.cc)
foreach(DIR ${FILES})
  get_filename_component(DIR "${DIR}" DIRECTORY)
  list(APPEND DIRS ${DIR})
endforeach(DIR)
list(REMOVE_DUPLICATES DIRS)
set(CSG_RST_FILES)
foreach(DIR ${DIRS})
  file(GLOB SRC ${DIR}/*.cc)
  add_executable(${DIR} ${SRC})
  target_link_libraries(${DIR} VOTCA::votca_csg)
  if (BUILD_MANPAGES)
    add_custom_command(OUTPUT ${DIR}.man
      COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_help2doc> --name ${DIR}/csg_${DIR} --format groff --out ${DIR}.man
      COMMENT "Building csg_${DIR} manpage"
      DEPENDS VOTCA::votca_help2doc ${DIR})
    add_custom_target(${DIR}_manpage DEPENDS ${DIR}.man)
    add_dependencies(manpages ${DIR}_manpage)
  endif()
  if(SPHINX_FOUND)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/csg/csg_${DIR}.rst
      COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_help2doc> --name ${DIR}/csg_${DIR} --format rst --out ${VOTCA_SPHINX_DIR}/csg/csg_${DIR}.rst
      COMMENT "Building csg_${DIR} rst doc"
      DEPENDS VOTCA::votca_help2doc ${DIR})
    list(APPEND CSG_RST_FILES ${VOTCA_SPHINX_DIR}/csg/csg_${DIR}.rst)
  endif()
  if(INSTALL_CSGAPPS)
    install(TARGETS ${DIR} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    if(BUILD_MANPAGES)
      install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${DIR}.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME csg_${DIR}.1)
    endif()
  endif()
  set_target_properties(${DIR} PROPERTIES OUTPUT_NAME csg_${DIR} RUNTIME_OUTPUT_DIRECTORY ${DIR})
  add_test(${DIR}Help ${DIR}/csg_${DIR} --help)
endforeach(DIR)

if(SPHINX_FOUND)
  add_custom_target(csgapps-rst DEPENDS ${CSG_RST_FILES})
  add_dependencies(sphinx-deps csgapps-rst)
endif()
