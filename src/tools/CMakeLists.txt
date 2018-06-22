foreach(PROG csg_reupdate csg_map csg_dump csg_property csg_resample csg_stat csg_fmatch csg_imcrepack csg_gmxtopol csg_dlptopol csg_density)
  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  target_link_libraries(${PROG} votca_csg)
  install(TARGETS ${PROG} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  if (TXT2TAGS_FOUND AND BASH)
    add_custom_command(OUTPUT ${PROG}.man
      COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PROG} --help > ${PROG}.help
      COMMAND ${BASH} ${CMAKE_CURRENT_BINARY_DIR}/../../scripts/help2t2t ${PROG}.help > ${PROG}.t2t
      COMMAND ${TXT2TAGS_EXECUTABLE} -q -t man -i ${PROG}.t2t -o ${PROG}.man
      DEPENDS help2t2t_build ${PROG})
    add_custom_target(${PROG}_manpage DEPENDS ${PROG}.man)
    add_dependencies(manpages ${PROG}_manpage)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROG}.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME ${PROG}.1)
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${PROG}.help ${PROG}.t2t)
  endif(TXT2TAGS_FOUND AND BASH)

  if(ENABLE_TESTING)
    add_test(${PROG}Help ${PROG} --help)
    # run tests for tools (for coverage) as well
    set_tests_properties(${PROG}Help PROPERTIES LABELS "csg;tools;votca")
  endif(ENABLE_TESTING)
endforeach(PROG)

if(ENABLE_TESTING)
  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_reupdate)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/spce)
  file(MAKE_DIRECTORY ${RUNPATH})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.param.in_re ${RUNPATH}/CG-CG.param.cur)
  add_test(NAME Run_csg_reupdate COMMAND csg_reupdate --gentable true --options ${REFPATH}/settings_re.xml WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_reupdate PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_reupdate_output COMMAND ${VOTCA_COMPARE} -f1 CG-CG.pot.new -f2 ${REFPATH}/CG-CG.pot.re WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_reupdate_output PROPERTIES DEPENDS Run_csg_reupdate)
  set_tests_properties(Compare_csg_reupdate_output PROPERTIES LABELS "csg;tools;votca")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_fmatch)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME Run_csg_fmatch 
    COMMAND csg_fmatch --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --options ${REFPATH}/settings_fmatch.xml
                      --cg ${REFPATH}/mapping.xml 
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_fmatch PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_fmatch_output COMMAND ${VOTCA_COMPARE} -f1 CG-CG.force -f2 ${REFPATH}/CG-CG.force.fmatch WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_fmatch_output PROPERTIES DEPENDS Run_csg_fmatch)
  set_tests_properties(Compare_csg_fmatch_output PROPERTIES LABELS "csg;tools;votca")
endif()
