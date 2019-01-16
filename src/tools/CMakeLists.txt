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
  add_test(NAME Compare_csg_reupdate_output COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.pot.new -f2 ${REFPATH}/CG-CG.pot.re WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_reupdate_output PROPERTIES DEPENDS Run_csg_reupdate)
  set_tests_properties(Compare_csg_reupdate_output PROPERTIES LABELS "csg;tools;votca")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_fmatch)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME Run_csg_fmatch 
    COMMAND csg_fmatch --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --options ${REFPATH}/settings_fmatch.xml
                      --cg ${REFPATH}/mapping.xml 
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_fmatch PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_fmatch_output COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.force -f2 ${REFPATH}/CG-CG.force.fmatch WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_fmatch_output PROPERTIES DEPENDS Run_csg_fmatch)
  set_tests_properties(Compare_csg_fmatch_output PROPERTIES LABELS "csg;tools;votca")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_fmatch_3body)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME Run_csg_fmatch_3body 
    COMMAND csg_fmatch --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --options ${REFPATH}/settings_fmatch_3body.xml
                      --cg ${REFPATH}/mapping.xml 
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_fmatch_3body PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_fmatch_3body_output1 COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.force -f2 ${REFPATH}/CG-CG.force.fmatch_3body WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME Compare_csg_fmatch_3body_output2 COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG-CG.force -f2 ${REFPATH}/CG-CG-CG.force.fmatch_3body WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME Compare_csg_fmatch_3body_output3 COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG-CG.pot -f2 ${REFPATH}/CG-CG-CG.pot.fmatch_3body WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_fmatch_3body_output1 PROPERTIES DEPENDS Run_csg_fmatch)
  set_tests_properties(Compare_csg_fmatch_3body_output1 PROPERTIES LABELS "csg;tools;votca")
  set_tests_properties(Compare_csg_fmatch_3body_output2 PROPERTIES DEPENDS Run_csg_fmatch)
  set_tests_properties(Compare_csg_fmatch_3body_output2 PROPERTIES LABELS "csg;tools;votca")
  set_tests_properties(Compare_csg_fmatch_3body_output3 PROPERTIES DEPENDS Run_csg_fmatch)
  set_tests_properties(Compare_csg_fmatch_3body_output3 PROPERTIES LABELS "csg;tools;votca")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_map)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME Run_csg_map
    COMMAND csg_map --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --out frame.dump --no-map --force
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_map PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_map_output COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 frame.dump -f2 ${REFPATH}/frame.dump WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_map_output PROPERTIES DEPENDS Run_csg_map)
  set_tests_properties(Compare_csg_map_output PROPERTIES LABELS "csg;tools;votca")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_stat)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME Run_csg_stat
    COMMAND csg_stat --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump
                     --options  ${REFPATH}/settings_rdf.xml --cg ${REFPATH}/mapping.xml
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_stat PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_stat_output1 COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.dist.new -f2 ${REFPATH}/CG-CG.rdf WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME Compare_csg_stat_output2 COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.force.new -f2 ${REFPATH}/CG-CG.pmf WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_stat_output1 PROPERTIES DEPENDS Run_csg_stat)
  set_tests_properties(Compare_csg_stat_output1 PROPERTIES LABELS "csg;tools;votca")
  set_tests_properties(Compare_csg_stat_output2 PROPERTIES DEPENDS Run_csg_stat)
  set_tests_properties(Compare_csg_stat_output2 PROPERTIES LABELS "csg;tools;votca") 

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_stat_angular)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME Run_csg_stat_angular
    COMMAND csg_stat --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump
                     --options  ${REFPATH}/settings_angular.xml --cg ${REFPATH}/mapping.xml
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_stat_angular PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_stat_angular_output COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG-CG.dist.new -f2 ${REFPATH}/CG-CG-CG.angular WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_stat_angular_output PROPERTIES DEPENDS Run_csg_stat_angular)
  set_tests_properties(Compare_csg_stat_angular_output PROPERTIES LABELS "csg;tools;votca")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_stat-imc)
  file(MAKE_DIRECTORY ${RUNPATH})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.imc.tgt ${RUNPATH}/CG-CG.dist.tgt)
  add_test(NAME Run_csg_stat-imc
    COMMAND csg_stat --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --do-imc
                     --options  ${REFPATH}/settings_imc.xml --cg ${REFPATH}/mapping.xml
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_stat-imc PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_stat-imc_output COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.imc -f2 ${REFPATH}/CG-CG.imc WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_stat-imc_output PROPERTIES DEPENDS Run_csg_stat-imc)
  set_tests_properties(Compare_csg_stat-imc_output PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_stat-imc_output_2 COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.dist.new -f2 ${REFPATH}/CG-CG.dist.imc WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_stat-imc_output_2 PROPERTIES DEPENDS Run_csg_stat-imc)
  set_tests_properties(Compare_csg_stat-imc_output_2 PROPERTIES LABELS "csg;tools;votca")
  
  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_resample)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/csg_resample)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME Run_csg_resample_akima
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_akima --type akima --grid 0:0.1:3
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_resample_akima PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_resample_akima_output COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 table_akima -f2 ${REFPATH}/table_akima WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_resample_akima_output PROPERTIES DEPENDS Run_csg_resample_akima)
  set_tests_properties(Compare_csg_resample_akima_output PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Run_csg_resample_cubic
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_cubic --type cubic --grid 0:0.1:3
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_resample_cubic PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_resample_cubic_output COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 table_cubic -f2 ${REFPATH}/table_cubic WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_resample_cubic_output PROPERTIES DEPENDS Run_csg_resample_cubic)
  set_tests_properties(Compare_csg_resample_cubic_output PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Run_csg_resample_linear
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_linear --type linear --grid 0:0.1:3
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_resample_linear PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_resample_linear_output COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 table_linear -f2 ${REFPATH}/table_linear WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_resample_linear_output PROPERTIES DEPENDS Run_csg_resample_linear)
  set_tests_properties(Compare_csg_resample_linear_output PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Run_csg_resample_cubicfit
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_cubicfit --type cubic --grid 0:1:9 --fitgrid 0:2:9
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_resample_cubicfit PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_resample_cubicfit_output COMMAND ${VOTCA_COMPARE} --etol ${REGRESSIONTEST_TOLERANCE} -f1 table_cubicfit -f2 ${REFPATH}/table_cubicfit WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_resample_cubicfit_output PROPERTIES DEPENDS Run_csg_resample_cubicfit)
  set_tests_properties(Compare_csg_resample_cubicfit_output PROPERTIES LABELS "csg;tools;votca")
  
  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_dump)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/methanol-water/)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME Run_csg_dump
    COMMAND ${BASH} -c "$<TARGET_FILE:csg_dump> --top ${REFPATH}/topol_cg.xml > csg_dump.out"
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Run_csg_dump PROPERTIES LABELS "csg;tools;votca")
  add_test(NAME Compare_csg_dump_output COMMAND ${CMAKE_COMMAND} -E compare_files csg_dump.out ${REFPATH}/csg_dump.out WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(Compare_csg_dump_output PROPERTIES DEPENDS Run_csg_dump)
  set_tests_properties(Compare_csg_dump_output PROPERTIES LABELS "csg;tools;votca")
endif()
