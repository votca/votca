set(CSG_RST_FILES)
foreach(PROG csg_reupdate csg_map csg_dump csg_property csg_resample csg_stat csg_fmatch csg_gmxtopol csg_dlptopol csg_density csg_imc_solve)
  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  target_link_libraries(${PROG} votca_csg)
  install(TARGETS ${PROG} EXPORT VOTCA_CSG_Targets RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  add_executable(VOTCA::${PROG} ALIAS ${PROG})
  set_property(TARGET votca_csg APPEND PROPERTY BINARIES "${PROG}")
  if (BUILD_MANPAGES)
    add_custom_command(OUTPUT ${PROG}.man
      COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_help2doc> --name $<TARGET_FILE:VOTCA::${PROG}> --format groff --out ${PROG}.man
      COMMENT "Building ${PROG} manpage"
      DEPENDS VOTCA::votca_help2doc ${PROG})
    add_custom_target(${PROG}_manpage DEPENDS ${PROG}.man)
    add_dependencies(manpages ${PROG}_manpage)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROG}.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME ${PROG}.1)
  endif()
  if(SPHINX_FOUND)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/csg/${PROG}.rst
      COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_help2doc> --name $<TARGET_FILE:VOTCA::${PROG}> --format rst --out ${VOTCA_SPHINX_DIR}/csg/${PROG}.rst
      COMMENT "Building ${PROG} rst doc"
      DEPENDS VOTCA::votca_help2doc ${PROG})
    list(APPEND CSG_RST_FILES ${VOTCA_SPHINX_DIR}/csg/${PROG}.rst)
  endif()

  add_test(integration_${PROG}Help ${PROG} --help)
endforeach(PROG)

if(SPHINX_FOUND)
  add_custom_target(csg-tools-rst DEPENDS ${CSG_RST_FILES})
  add_dependencies(sphinx-deps csg-tools-rst)
endif()

if(BUILD_TESTING)
  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_reupdate)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/spce)
  file(MAKE_DIRECTORY ${RUNPATH})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.param.in_re ${RUNPATH}/CG-CG.param.cur)
  add_test(NAME integration_Run_csg_reupdate COMMAND csg_reupdate --gentable true --options ${REFPATH}/settings_re.xml WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_reupdate_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.pot.new -f2 ${REFPATH}/CG-CG.pot.re WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_reupdate_output PROPERTIES DEPENDS integration_Run_csg_reupdate)

  
  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_reupdate2)
  file(MAKE_DIRECTORY ${RUNPATH})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.rdf ${RUNPATH}/CG-CG.dist.new)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.imc.tgt ${RUNPATH}/CG-CG.dist.tgt)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.param.in_re ${RUNPATH}/CG-CG.param.cur)
  add_test(NAME integration_Run_csg_reupdate2 COMMAND csg_reupdate --options ${REFPATH}/settings_re.xml --top ${REFPATH}/topol_cg.xml --trj ${REFPATH}/frame_cg.dump --hessian-check no WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_reupdate_output2 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.pot.new -f2 ${REFPATH}/CG-CG.pot.new.re WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_reupdate_output3 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.param.new -f2 ${REFPATH}/CG-CG.param.re WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_reupdate_output2 PROPERTIES DEPENDS integration_Run_csg_reupdate2)
  set_tests_properties(integration_Compare_csg_reupdate_output3 PROPERTIES DEPENDS integration_Run_csg_reupdate2)

  if(GMX_FOUND AND GROMACS_FOUND)
    set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_reupdate_gmx)
    file(MAKE_DIRECTORY ${RUNPATH})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/topol_cg.top ${RUNPATH}/topol.top)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/grompp_cg.mdp ${RUNPATH}/grompp.mdp)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/conf_cg.gro ${RUNPATH}/conf.gro)
    add_test(NAME integration_Run_gmx_grompp_re COMMAND ${GMX_EXECUTABLE} grompp WORKING_DIRECTORY ${RUNPATH})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.rdf ${RUNPATH}/CG-CG.dist.new)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.imc.tgt ${RUNPATH}/CG-CG.dist.tgt)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.param.in_re ${RUNPATH}/CG-CG.param.cur)
    add_test(NAME integration_Run_csg_reupdate_gmx COMMAND csg_reupdate --options ${REFPATH}/settings_re.xml --top topol.tpr --trj conf.gro --hessian-check no WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Run_csg_reupdate_gmx PROPERTIES DEPENDS integration_Run_gmx_grompp_re)
    add_test(NAME integration_Compare_csg_reupdate_gmx_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.pot.new -f2 ${REFPATH}/CG-CG.pot.new.re_gmx WORKING_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Compare_csg_reupdate_gmx_output2 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.param.new -f2 ${REFPATH}/CG-CG.param.re_gmx WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Compare_csg_reupdate_gmx_output PROPERTIES DEPENDS integration_Run_csg_reupdate_gmx)
    set_tests_properties(integration_Compare_csg_reupdate_gmx_output2 PROPERTIES DEPENDS integration_Run_csg_reupdate_gmx)
  endif()

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_fmatch)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_fmatch 
    COMMAND csg_fmatch --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --options ${REFPATH}/settings_fmatch.xml
                      --cg ${REFPATH}/mapping.xml 
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_fmatch_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.force -f2 ${REFPATH}/CG-CG.force.fmatch WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_fmatch_output PROPERTIES DEPENDS integration_Run_csg_fmatch)

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_fmatch_3body)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_fmatch_3body 
    COMMAND csg_fmatch --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --options ${REFPATH}/settings_fmatch_3body.xml
                      --cg ${REFPATH}/mapping.xml 
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_fmatch_3body_output1 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.force -f2 ${REFPATH}/CG-CG.force.fmatch_3body WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_fmatch_3body_output2 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG-CG.force -f2 ${REFPATH}/CG-CG-CG.force.fmatch_3body WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_fmatch_3body_output3 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG-CG.pot -f2 ${REFPATH}/CG-CG-CG.pot.fmatch_3body WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_fmatch_3body_output1 PROPERTIES DEPENDS integration_Run_csg_fmatch_3body)
  set_tests_properties(integration_Compare_csg_fmatch_3body_output2 PROPERTIES DEPENDS integration_Run_csg_fmatch_3body)
  set_tests_properties(integration_Compare_csg_fmatch_3body_output3 PROPERTIES DEPENDS integration_Run_csg_fmatch_3body)

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_map)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_map
    COMMAND csg_map --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --out frame.dump --no-map --force
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_map_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 frame.dump -f2 ${REFPATH}/frame.dump WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_map_output PROPERTIES DEPENDS integration_Run_csg_map)

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_map_gro)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_map_gro COMMAND csg_map --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump  --cg ${REFPATH}/mapping.xml --out conf_cg.gro --nframes 1 WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_map_gro_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 conf_cg.gro -f2 ${REFPATH}/conf_cg.gro WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_map_gro_output PROPERTIES DEPENDS integration_Run_csg_map_gro)

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_stat)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_stat
    COMMAND csg_stat --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump
                     --options  ${REFPATH}/settings_rdf.xml --cg ${REFPATH}/mapping.xml
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_stat_output1 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.dist.new -f2 ${REFPATH}/CG-CG.rdf WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_stat_output2 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.force.new -f2 ${REFPATH}/CG-CG.pmf WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat_output1 PROPERTIES DEPENDS integration_Run_csg_stat)
  set_tests_properties(integration_Compare_csg_stat_output2 PROPERTIES DEPENDS integration_Run_csg_stat)

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_stat_angular)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_stat_angular
    COMMAND csg_stat --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump
                     --options  ${REFPATH}/settings_angular.xml --cg ${REFPATH}/mapping.xml
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_stat_angular_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG-CG.dist.new -f2 ${REFPATH}/CG-CG-CG.angular WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat_angular_output PROPERTIES DEPENDS integration_Run_csg_stat_angular)

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_stat-imc)
  file(MAKE_DIRECTORY ${RUNPATH})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/CG-CG.imc.tgt ${RUNPATH}/CG-CG.dist.tgt)
  add_test(NAME integration_Run_csg_stat-imc
    COMMAND csg_stat --top ${REFPATH}/topol.xml --trj ${REFPATH}/frame.dump --do-imc
                     --options  ${REFPATH}/settings_imc.xml --cg ${REFPATH}/mapping.xml
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_stat-imc_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.imc -f2 ${REFPATH}/CG-CG.imc WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat-imc_output PROPERTIES DEPENDS integration_Run_csg_stat-imc)
  add_test(NAME integration_Compare_csg_stat-imc_output_2 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.dist.new -f2 ${REFPATH}/CG-CG.dist.imc WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat-imc_output_2 PROPERTIES DEPENDS integration_Run_csg_stat-imc)
  add_test(NAME integration_Compare_csg_stat-imc_output_3 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.gmc -f2 ${REFPATH}/CG-CG.gmc WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat-imc_output_3 PROPERTIES DEPENDS integration_Run_csg_stat-imc)
  add_test(NAME integration_Compare_csg_stat-imc_output_4 COMMAND ${CMAKE_COMMAND} -E compare_files CG-CG.idx ${REFPATH}/CG-CG.idx WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat-imc_output_4 PROPERTIES DEPENDS integration_Run_csg_stat-imc)

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_imc_solve)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/spce/imc)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_imc_solve
    COMMAND csg_imc_solve -i ${REFPATH}/CG-CG.imc -g ${REFPATH}/CG-CG.gmc -n ${REFPATH}/CG-CG.idx
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_imc_solve_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 CG-CG.dpot.imc -f2 ${REFPATH}/CG-CG.dpot.imc
 WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_imc_solve_output PROPERTIES DEPENDS integration_Run_csg_imc_solve)

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_imc_solve_reg)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/LJ1-LJ2)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_imc_solve_reg
    COMMAND csg_imc_solve -r 1000 -i ${REFPATH}/group_1.imc -g ${REFPATH}/group_1.gmc -n ${REFPATH}/group_1.idx
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_imc_solve_reg_output_2 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> 
    --etol ${INTEGRATIONTEST_TOLERANCE} -f1 LJ1-LJ1.dpot.imc -f2 ${REFPATH}/LJ1-LJ1.dpot.imc
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_imc_solve_reg_output_2 PROPERTIES DEPENDS integration_Run_csg_imc_solve_reg)
  add_test(NAME integration_Compare_csg_imc_solve_reg_output_3 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> 
    --etol ${INTEGRATIONTEST_TOLERANCE} -f1 LJ1-LJ2.dpot.imc -f2 ${REFPATH}/LJ1-LJ2.dpot.imc
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_imc_solve_reg_output_3 PROPERTIES DEPENDS integration_Run_csg_imc_solve_reg)
  add_test(NAME integration_Compare_csg_imc_solve_reg_output_4 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> 
    --etol ${INTEGRATIONTEST_TOLERANCE} -f1 LJ2-LJ2.dpot.imc -f2 ${REFPATH}/LJ2-LJ2.dpot.imc
    WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_imc_solve_reg_output_4 PROPERTIES DEPENDS integration_Run_csg_imc_solve_reg)

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_stat_imc_multi)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/LJ1-LJ2/smaller_system)
  file(MAKE_DIRECTORY ${RUNPATH})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/LJ1-LJ1.dist.tgt ${RUNPATH}/LJ1-LJ1.dist.tgt)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/LJ1-LJ2.dist.tgt ${RUNPATH}/LJ1-LJ2.dist.tgt)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/LJ2-LJ2.dist.tgt ${RUNPATH}/LJ2-LJ2.dist.tgt)
  add_test(NAME integration_Run_csg_stat_imc_multi
    COMMAND csg_stat --top ${REFPATH}/topol.xml --trj ${REFPATH}/traj.gro --do-imc
                   --options ${REFPATH}/settings_imc.xml
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_stat-imc_multi_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 all.imc -f2 ${REFPATH}/all.imc WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat-imc_multi_output PROPERTIES DEPENDS integration_Run_csg_stat_imc_multi)
  add_test(NAME integration_Compare_csg_stat-imc_multi_output_2 COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 all.gmc -f2 ${REFPATH}/all.gmc WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat-imc_multi_output_2 PROPERTIES DEPENDS integration_Run_csg_stat_imc_multi)
  add_test(NAME integration_Compare_csg_stat-imc_multi_output_3 COMMAND ${CMAKE_COMMAND} -E compare_files all.idx ${REFPATH}/all.idx WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_stat-imc_multi_output_3 PROPERTIES DEPENDS integration_Run_csg_stat_imc_multi)

  set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_resample)
  set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/csg_resample)
  file(MAKE_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Run_csg_resample_akima
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_akima --type akima --grid 0:0.1:3
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_resample_akima_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 table_akima -f2 ${REFPATH}/table_akima WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_resample_akima_output PROPERTIES DEPENDS integration_Run_csg_resample_akima)
  add_test(NAME integration_Run_csg_resample_cubic
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_cubic --type cubic --grid 0:0.1:3
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_resample_cubic_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 table_cubic -f2 ${REFPATH}/table_cubic WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_resample_cubic_output PROPERTIES DEPENDS integration_Run_csg_resample_cubic)
  add_test(NAME integration_Run_csg_resample_linear
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_linear --type linear --grid 0:0.1:3
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_resample_linear_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 table_linear -f2 ${REFPATH}/table_linear WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_resample_linear_output PROPERTIES DEPENDS integration_Run_csg_resample_linear)
  add_test(NAME integration_Run_csg_resample_cubicfit
    COMMAND csg_resample --in ${REFPATH}/table_in --out table_cubicfit --type cubic --grid 0:1:9 --fitgrid 0:2:9
    WORKING_DIRECTORY ${RUNPATH})
  add_test(NAME integration_Compare_csg_resample_cubicfit_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 table_cubicfit -f2 ${REFPATH}/table_cubicfit WORKING_DIRECTORY ${RUNPATH})
  set_tests_properties(integration_Compare_csg_resample_cubicfit_output PROPERTIES DEPENDS integration_Run_csg_resample_cubicfit)
 
  if(BASH)
    set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_dump)
    set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/methanol-water/)
    file(MAKE_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Run_csg_dump
      COMMAND ${BASH} -c "$<TARGET_FILE:csg_dump> --top ${REFPATH}/topol_cg.xml > csg_dump.out"
      WORKING_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Compare_csg_dump_output COMMAND ${CMAKE_COMMAND} -E compare_files csg_dump.out ${REFPATH}/csg_dump.out WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Compare_csg_dump_output PROPERTIES DEPENDS integration_Run_csg_dump)

    set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_dump_gro)
    set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/spce/3mols)
    file(MAKE_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Run_csg_dump_gro
      COMMAND ${BASH} -c "$<TARGET_FILE:csg_dump> --top ${REFPATH}/conf.gro > csg_dump_gro.out"
      WORKING_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Compare_csg_dump_gro_output COMMAND ${CMAKE_COMMAND} -E compare_files csg_dump_gro.out ${REFPATH}/csg_dump_gro.out WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Compare_csg_dump_gro_output PROPERTIES DEPENDS integration_Run_csg_dump_gro)
    
    set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_dump_data)
    set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/spce/3mols_cg)
    file(MAKE_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Run_csg_dump_data
      COMMAND ${BASH} -c "$<TARGET_FILE:csg_dump> --top ${REFPATH}/spce.data > csg_dump_data.out"
      WORKING_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Compare_csg_dump_data_output COMMAND ${CMAKE_COMMAND} -E compare_files csg_dump_data.out ${REFPATH}/csg_dump_data.out WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Compare_csg_dump_data_output PROPERTIES DEPENDS integration_Run_csg_dump_data)
  endif()

  if(GMX_FOUND AND GROMACS_FOUND)
    set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_gmx_3mols)
    set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/spce/3mols)
    file(MAKE_DIRECTORY ${RUNPATH})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/topol.top ${RUNPATH}/topol.top)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/grompp.mdp ${RUNPATH}/grompp.mdp)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/conf.gro ${RUNPATH}/conf.gro)
    add_test(NAME integration_Run_gmx_grompp COMMAND ${GMX_EXECUTABLE} grompp WORKING_DIRECTORY ${RUNPATH})
    if(BASH)
      add_test(NAME integration_Run_csg_dump_gmx_tpr
        COMMAND ${BASH} -c "$<TARGET_FILE:csg_dump> --top topol.tpr > csg_dump_tpr.out"
        WORKING_DIRECTORY ${RUNPATH})
      set_tests_properties(integration_Run_csg_dump_gmx_tpr PROPERTIES DEPENDS integration_Run_gmx_grompp)
      add_test(NAME integration_Compare_csg_dump_gmx_tpr_output COMMAND ${CMAKE_COMMAND} -E compare_files csg_dump_tpr.out ${REFPATH}/csg_dump_tpr.out WORKING_DIRECTORY ${RUNPATH})
      set_tests_properties(integration_Compare_csg_dump_gmx_tpr_output PROPERTIES DEPENDS integration_Run_csg_dump_gmx_tpr)
      add_test(NAME integration_Run_csg_dump_gmx_tpr_excl
        COMMAND ${BASH} -c "$<TARGET_FILE:csg_dump> --excl --top topol.tpr > csg_dump_tpr_excl.out"
        WORKING_DIRECTORY ${RUNPATH})
      set_tests_properties(integration_Run_csg_dump_gmx_tpr_excl PROPERTIES DEPENDS integration_Run_gmx_grompp)
      add_test(NAME integration_Compare_csg_dump_gmx_tpr_excl_output COMMAND ${CMAKE_COMMAND} -E compare_files csg_dump_tpr_excl.out ${REFPATH}/csg_dump_tpr_excl.out WORKING_DIRECTORY ${RUNPATH})
      set_tests_properties(integration_Compare_csg_dump_gmx_tpr_excl_output PROPERTIES DEPENDS integration_Run_csg_dump_gmx_tpr_excl)
      add_test(NAME integration_Run_csg_dump_gmx_tpr_cg
        COMMAND ${BASH} -c "$<TARGET_FILE:csg_dump> --cg ${REFPATH}/../mapping.xml --top topol.tpr > csg_dump_tpr_cg.out"
        WORKING_DIRECTORY ${RUNPATH})
      set_tests_properties(integration_Run_csg_dump_gmx_tpr_cg PROPERTIES DEPENDS integration_Run_gmx_grompp)
      add_test(NAME integration_Compare_csg_dump_gmx_tpr_cg_output COMMAND ${CMAKE_COMMAND} -E compare_files csg_dump_tpr_cg.out ${REFPATH}/csg_dump_tpr_cg.out WORKING_DIRECTORY ${RUNPATH})
      set_tests_properties(integration_Compare_csg_dump_gmx_tpr_cg_output PROPERTIES DEPENDS integration_Run_csg_dump_gmx_tpr_cg)
    endif()
    add_test(NAME integration_Run_gmx_mdrun COMMAND ${GMX_EXECUTABLE} mdrun WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Run_gmx_mdrun PROPERTIES DEPENDS integration_Run_gmx_grompp)
    add_test(NAME integration_Run_csg_map_gmx COMMAND csg_map --top topol.tpr --trj traj.trr --out traj.gro --no-map --nframes 1 WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Run_csg_map_gmx PROPERTIES DEPENDS integration_Run_gmx_mdrun)
    add_test(NAME integration_Compare_csg_map_gmx_output COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 traj.gro -f2 ${REFPATH}/conf.gro WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Compare_csg_map_gmx_output PROPERTIES DEPENDS integration_Run_csg_map_gmx)
    add_test(NAME integration_Run_csg_map_gmx_xtc COMMAND csg_map --top topol.tpr --trj traj_comp.xtc  --out traj_comp.gro --no-map --nframes 1 WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Run_csg_map_gmx_xtc PROPERTIES DEPENDS integration_Run_gmx_mdrun)
    add_test(NAME integration_Compare_csg_map_gmx_xtc_output
             COMMAND ${Python_EXECUTABLE} $<TARGET_FILE:VOTCA::votca_compare> --etol ${INTEGRATIONTEST_TOLERANCE} -f1 traj_comp.gro -f2 ${REFPATH}/conf.gro WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Compare_csg_map_gmx_xtc_output PROPERTIES DEPENDS integration_Run_csg_map_gmx_xtc)
  endif()

  if(BASH)
    set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_property)
    set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/csg_property)
    file(MAKE_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Run_csg_property
      COMMAND ${BASH} -c "$<TARGET_FILE:csg_property> --file ${REFPATH}/csg_defaults.xml --print . --short --path cg.inverse.gromacs.DESC > csg_property.out"
      WORKING_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Compare_csg_property_output COMMAND ${CMAKE_COMMAND} -E compare_files csg_property.out ${REFPATH}/csg_property.out WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Compare_csg_property_output PROPERTIES DEPENDS integration_Run_csg_property)
    add_test(NAME integration_Run_csg_property2
      COMMAND ${BASH} -c "$<TARGET_FILE:csg_property> --file ${REFPATH}/csg_defaults.xml --print . --path cg.inverse.gromacs.DESC > csg_property.out2"
      WORKING_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Compare_csg_property_output2 COMMAND ${CMAKE_COMMAND} -E compare_files csg_property.out2 ${REFPATH}/csg_property.out2 WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Compare_csg_property_output2 PROPERTIES DEPENDS integration_Run_csg_property2)
    add_test(NAME integration_Run_csg_property3
      COMMAND ${BASH} -c "$<TARGET_FILE:csg_property> --file ${REFPATH}/test.xml --path A.B.C --print . > csg_property.out3"
      WORKING_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Compare_csg_property_output3 COMMAND ${CMAKE_COMMAND} -E compare_files csg_property.out3 ${REFPATH}/csg_property.out3 WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Compare_csg_property_output3 PROPERTIES DEPENDS integration_Run_csg_property3)
    add_test(NAME integration_Run_csg_property4
      COMMAND ${BASH} -c "$<TARGET_FILE:csg_property> --file ${REFPATH}/test.xml --path A.B --print C > csg_property.out4"
      WORKING_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Compare_csg_property_output4 COMMAND ${CMAKE_COMMAND} -E compare_files csg_property.out4 ${REFPATH}/csg_property.out4 WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Compare_csg_property_output4 PROPERTIES DEPENDS integration_Run_csg_property4)

  endif(BASH)

  if(GMX_FOUND AND GROMACS_FOUND)
    set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_boltzmann_ppy)
    set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/ppy/)
    file(MAKE_DIRECTORY ${RUNPATH})
    foreach(_SRC topol.top grompp.mdp conf.gro ffgmx.itp ppy10.xml)
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/${_SRC} ${RUNPATH}/${_SRC})
    endforeach()
    add_test(NAME integration_Run_gmx_grompp_ppy COMMAND ${GMX_EXECUTABLE} grompp WORKING_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Run_csg_boltzmann_ppy
      COMMAND $<TARGET_FILE:csg_boltzmann> --top topol.tpr --cg ppy10.xml --excl excl.txt
      WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Run_csg_boltzmann_ppy PROPERTIES DEPENDS integration_Run_gmx_grompp_ppy)
    add_test(NAME Compare_csg_boltzmann_output_ppy COMMAND ${CMAKE_COMMAND} -E compare_files excl.txt ${REFPATH}/excl.txt WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(Compare_csg_boltzmann_output_ppy PROPERTIES DEPENDS integration_Run_csg_boltzmann_ppy)
    
    set(RUNPATH ${CMAKE_CURRENT_BINARY_DIR}/Run_csg_boltzmann_mma)
    set(REFPATH ${CMAKE_CURRENT_SOURCE_DIR}/references/mma/)
    file(MAKE_DIRECTORY ${RUNPATH})
    foreach(_SRC topol.top grompp.mdp conf.gro ffgmx.itp mapping.xml)
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${REFPATH}/${_SRC} ${RUNPATH}/${_SRC})
    endforeach()
    add_test(NAME integration_Run_gmx_grompp_mma COMMAND ${GMX_EXECUTABLE} grompp WORKING_DIRECTORY ${RUNPATH})
    add_test(NAME integration_Run_csg_boltzmann_mma
      COMMAND $<TARGET_FILE:csg_boltzmann> --top topol.tpr --cg mapping.xml --excl excl.txt
      WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(integration_Run_csg_boltzmann_mma PROPERTIES DEPENDS integration_Run_gmx_grompp_mma)
    add_test(NAME Compare_csg_boltzmann_output_mma COMMAND ${CMAKE_COMMAND} -E compare_files excl.txt ${REFPATH}/excl.txt WORKING_DIRECTORY ${RUNPATH})
    set_tests_properties(Compare_csg_boltzmann_output_mma PROPERTIES DEPENDS integration_Run_csg_boltzmann_mma)
  endif(GMX_FOUND AND GROMACS_FOUND)
endif()
