foreach(PROG csg_reupdate csg_map csg_dump csg_property csg_resample csg_stat csg_fmatch csg_gmxtopol csg_dlptopol csg_density csg_imc_solve)
  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  target_link_libraries(${PROG} votca_csg)
  install(TARGETS ${PROG} EXPORT VOTCA_CSG_Targets RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(VOTCA::${PROG} ALIAS ${PROG})
  set_property(TARGET votca_csg APPEND PROPERTY BINARIES "${PROG}")
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
    add_test(integration_${PROG}Help ${PROG} --help)
    # run tests for tools (for coverage) as well
    set_tests_properties(integration_${PROG}Help PROPERTIES LABELS "csg;tools;votca;integration")
  endif(ENABLE_TESTING)
endforeach(PROG)

if(ENABLE_TESTING)
  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_reupdate)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/spce)
  file(MAKE_DIRECTORY ${RUNPATH})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.param.in_re ${RUNPATH}/CG-CG.param.cur)
  add_test(NAME integration_Run_csg_reupdate COMMAND csg_reupdate --gentable true --options ${REFPATH}/settings_re.xml WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_reupdate PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_reupdate_output COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.pot.new -f2 ${REFPATH}/CG-CG.pot.re WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_reupdate_output PROPERTIES DEPENDS integration_Run_csg_reupdate)
  set_tests_properties(integration_Compare_csg_reupdate_output PROPERTIES LABELS "csg;tools;votca;integration")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_fmatch)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_fmatch 
    COMMAND csg_fmatch --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --options ${REFPATH}/settings_fmatch.xml
                      --cg ${REFPATH}/mapping.xml 
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_fmatch PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_fmatch_output COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.force -f2 ${REFPATH}/CG-CG.force.fmatch WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_fmatch_output PROPERTIES DEPENDS integration_Run_csg_fmatch)
  set_tests_properties(integration_Compare_csg_fmatch_output PROPERTIES LABELS "csg;tools;votca;integration")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_fmatch_3body)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_fmatch_3body 
    COMMAND csg_fmatch --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --options ${REFPATH}/settings_fmatch_3body.xml
                      --cg ${REFPATH}/mapping.xml 
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_fmatch_3body PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_fmatch_3body_output1 COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.force -f2 ${REFPATH}/CG-CG.force.fmatch_3body WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_fmatch_3body_output2 COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG-CG.force -f2 ${REFPATH}/CG-CG-CG.force.fmatch_3body WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_fmatch_3body_output3 COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG-CG.pot -f2 ${REFPATH}/CG-CG-CG.pot.fmatch_3body WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_fmatch_3body_output1 PROPERTIES DEPENDS integration_Run_csg_fmatch)
  set_tests_properties(integration_Compare_csg_fmatch_3body_output1 PROPERTIES LABELS "csg;tools;votca;integration")
  set_tests_properties(integration_Compare_csg_fmatch_3body_output2 PROPERTIES DEPENDS integration_Run_csg_fmatch)
  set_tests_properties(integration_Compare_csg_fmatch_3body_output2 PROPERTIES LABELS "csg;tools;votca;integration")
  set_tests_properties(integration_Compare_csg_fmatch_3body_output3 PROPERTIES DEPENDS integration_Run_csg_fmatch)
  set_tests_properties(integration_Compare_csg_fmatch_3body_output3 PROPERTIES LABELS "csg;tools;votca;integration")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_map)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_map
    COMMAND csg_map --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --out frame.dump --no-map --force
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_map PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_map_output COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 frame.dump -f2 ${REFPATH}/frame.dump WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_map_output PROPERTIES DEPENDS integration_Run_csg_map)
  set_tests_properties(integration_Compare_csg_map_output PROPERTIES LABELS "csg;tools;votca;integration")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_stat)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_stat
    COMMAND csg_stat --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump
                     --options  ${REFPATH}/settings_rdf.xml --cg ${REFPATH}/mapping.xml
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_stat PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_stat_output1 COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.dist.new -f2 ${REFPATH}/CG-CG.rdf WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_stat_output2 COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.force.new -f2 ${REFPATH}/CG-CG.pmf WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat_output1 PROPERTIES DEPENDS integration_Run_csg_stat)
  set_tests_properties(integration_Compare_csg_stat_output1 PROPERTIES LABELS "csg;tools;votca;integration")
  set_tests_properties(integration_Compare_csg_stat_output2 PROPERTIES DEPENDS integration_Run_csg_stat)
  set_tests_properties(integration_Compare_csg_stat_output2 PROPERTIES LABELS "csg;tools;votca;integration") 

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_stat_angular)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_stat_angular
    COMMAND csg_stat --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump
                     --options  ${REFPATH}/settings_angular.xml --cg ${REFPATH}/mapping.xml
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_stat_angular PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_stat_angular_output COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG-CG.dist.new -f2 ${REFPATH}/CG-CG-CG.angular WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat_angular_output PROPERTIES DEPENDS integration_Run_csg_stat_angular)
  set_tests_properties(integration_Compare_csg_stat_angular_output PROPERTIES LABELS "csg;tools;votca;integration")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_stat-imc)
  file(MAKE_DIRECTORY ${RUNPATH})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.imc.tgt ${RUNPATH}/CG-CG.dist.tgt)
  add_test(NAME integration_Run_csg_stat-imc
    COMMAND csg_stat --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --do-imc
                     --options  ${REFPATH}/settings_imc.xml --cg ${REFPATH}/mapping.xml
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_stat-imc PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_stat-imc_output COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.imc -f2 ${REFPATH}/CG-CG.imc WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat-imc_output PROPERTIES DEPENDS integration_Run_csg_stat-imc)
  set_tests_properties(integration_Compare_csg_stat-imc_output PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_stat-imc_output_2 COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.dist.new -f2 ${REFPATH}/CG-CG.dist.imc WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat-imc_output_2 PROPERTIES DEPENDS integration_Run_csg_stat-imc)
  set_tests_properties(integration_Compare_csg_stat-imc_output_2 PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_stat-imc_output_3 COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 CG-CG.gmc -f2 ${REFPATH}/CG-CG.gmc WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat-imc_output_3 PROPERTIES DEPENDS integration_Run_csg_stat-imc)
  set_tests_properties(integration_Compare_csg_stat-imc_output_3 PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_stat-imc_output_4 COMMAND ${CMAKE_COMMAND} -E compare_files CG-CG.idx ${REFPATH}/CG-CG.idx WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat-imc_output_4 PROPERTIES DEPENDS integration_Run_csg_stat-imc)
  set_tests_properties(integration_Compare_csg_stat-imc_output_4 PROPERTIES LABELS "csg;tools;votca;integration")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_imc_solve)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/spce/imc)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_imc_solve
    COMMAND csg_imc_solve -i ${REFPATH}/CG-CG.imc -g ${REFPATH}/CG-CG.gmc -o table_csg_imc_solve.dat
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_imc_solve PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_imc_solve_output COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 table_csg_imc_solve.dat -f2 ${REFPATH}/CG-CG.sol
 WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_imc_solve_output PROPERTIES DEPENDS integration_Run_csg_imc_solve)
  set_tests_properties(integration_Compare_csg_imc_solve_output PROPERTIES LABELS "csg;tools;votca;integration")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_imc_solve_reg)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/LJ1-LJ2)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_imc_solve_reg
    COMMAND csg_imc_solve -r 1000 -i ${REFPATH}/group_1.imc -g ${REFPATH}/group_1.gmc
                          -n ${REFPATH}/group_1.idx -o group_1.sol
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_imc_solve_reg PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_imc_solve_reg_output COMMAND $<TARGET_FILE:VOTCA::votca_compare> 
    --etol ${REGRESSIONTEST_TOLERANCE} -f1 group_1.sol -f2 ${REFPATH}/group_1.sol
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_imc_solve_reg_output PROPERTIES DEPENDS integration_Run_csg_imc_solve_reg)
  set_tests_properties(integration_Compare_csg_imc_solve_reg_output PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_imc_solve_reg_output_2 COMMAND $<TARGET_FILE:VOTCA::votca_compare> 
    --etol ${REGRESSIONTEST_TOLERANCE} -f1 LJ1-LJ1.dpot.imc -f2 ${REFPATH}/LJ1-LJ1.dpot.imc
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_imc_solve_reg_output_2 PROPERTIES DEPENDS integration_Run_csg_imc_solve_reg)
  set_tests_properties(integration_Compare_csg_imc_solve_reg_output_2 PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_imc_solve_reg_output_3 COMMAND $<TARGET_FILE:VOTCA::votca_compare> 
    --etol ${REGRESSIONTEST_TOLERANCE} -f1 LJ1-LJ2.dpot.imc -f2 ${REFPATH}/LJ1-LJ2.dpot.imc
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_imc_solve_reg_output_3 PROPERTIES DEPENDS integration_Run_csg_imc_solve_reg)
  set_tests_properties(integration_Compare_csg_imc_solve_reg_output_3 PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_imc_solve_reg_output_4 COMMAND $<TARGET_FILE:VOTCA::votca_compare> 
    --etol ${REGRESSIONTEST_TOLERANCE} -f1 LJ2-LJ2.dpot.imc -f2 ${REFPATH}/LJ2-LJ2.dpot.imc
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_imc_solve_reg_output_4 PROPERTIES DEPENDS integration_Run_csg_imc_solve_reg)
  set_tests_properties(integration_Compare_csg_imc_solve_reg_output_4 PROPERTIES LABELS "csg;tools;votca;integration")

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_resample)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/csg_resample)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_resample_akima
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_akima --type akima --grid 0:0.1:3
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_resample_akima PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_resample_akima_output COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 table_akima -f2 ${REFPATH}/table_akima WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_resample_akima_output PROPERTIES DEPENDS integration_Run_csg_resample_akima)
  set_tests_properties(integration_Compare_csg_resample_akima_output PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Run_csg_resample_cubic
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_cubic --type cubic --grid 0:0.1:3
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_resample_cubic PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_resample_cubic_output COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 table_cubic -f2 ${REFPATH}/table_cubic WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_resample_cubic_output PROPERTIES DEPENDS integration_Run_csg_resample_cubic)
  set_tests_properties(integration_Compare_csg_resample_cubic_output PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Run_csg_resample_linear
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_linear --type linear --grid 0:0.1:3
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_resample_linear PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_resample_linear_output COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 table_linear -f2 ${REFPATH}/table_linear WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_resample_linear_output PROPERTIES DEPENDS integration_Run_csg_resample_linear)
  set_tests_properties(integration_Compare_csg_resample_linear_output PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Run_csg_resample_cubicfit
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_cubicfit --type cubic --grid 0:1:9 --fitgrid 0:2:9
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Run_csg_resample_cubicfit PROPERTIES LABELS "csg;tools;votca;integration")
  add_test(NAME integration_Compare_csg_resample_cubicfit_output COMMAND $<TARGET_FILE:VOTCA::votca_compare> --etol ${REGRESSIONTEST_TOLERANCE} -f1 table_cubicfit -f2 ${REFPATH}/table_cubicfit WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_resample_cubicfit_output PROPERTIES DEPENDS integration_Run_csg_resample_cubicfit)
  set_tests_properties(integration_Compare_csg_resample_cubicfit_output PROPERTIES LABELS "csg;tools;votca;integration")
 
  if(BASH)
    set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_dump)
    set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/methanol-water/)
    file(MAKE_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Run_csg_dump
      COMMAND ${BASH} -c "$<TARGET_FILE:csg_dump> --top ${REFPATH}/topol_cg.xml > csg_dump.out"
      WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Run_csg_dump PROPERTIES LABELS "csg;tools;votca;integration")
    add_test(NAME integration_Compare_csg_dump_output COMMAND ${CMAKE_COMMAND} -E compare_files csg_dump.out ${REFPATH}/csg_dump.out WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Compare_csg_dump_output PROPERTIES DEPENDS integration_Run_csg_dump)
    set_tests_properties(integration_Compare_csg_dump_output PROPERTIES LABELS "csg;tools;votca;integration")
  endif(BASH)
endif()
