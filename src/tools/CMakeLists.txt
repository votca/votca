foreach(PROG csg_map csg_dump csg_property csg_resample csg_stat csg_fmatch csg_imcrepack csg_gmxtopol csg_part_dist csg_density)
  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  target_link_libraries(${PROG} votca_csg)
  install(TARGETS ${PROG} RUNTIME DESTINATION bin)
  if (TXT2TAGS_FOUND AND BASH)
    add_custom_command(OUTPUT ${PROG}.man 
      COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PROG} --help > ${PROG}.help
      COMMAND ${BASH} ${CMAKE_BINARY_DIR}/scripts/help2t2t ${PROG}.help > ${PROG}.t2t 
      COMMAND ${TXT2TAGS_EXECUTABLE} -q -t man -i ${PROG}.t2t -o ${PROG}.man
      DEPENDS ${CMAKE_BINARY_DIR}/scripts/help2t2t ${PROG})
    add_custom_target(${PROG}_manpage DEPENDS ${PROG}.man)
    add_dependencies(manpages ${PROG}_manpage)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROG}.man DESTINATION ${MAN}/man1 RENAME ${PROG}.1)
  endif(TXT2TAGS_FOUND AND BASH)
endforeach(PROG)
