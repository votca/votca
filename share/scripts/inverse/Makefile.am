inversedir = ${pkgdatadir}/scripts/inverse

inverse_SCRIPTS = \
  add_pot_ibm.sh            functions_gromacs.sh  prepare_gromacs.sh   table_to_xvg.pl \
  add_POT.pl                functions_perl.pl     post_add.sh              pressure_cor.pl       \
  calc_pressure_gromacs.sh  initalize_gromacs.sh  post_add_single.sh       RDF_to_POT.pl        update_ibm_pot.pl \
  calc_rdf_gromacs.sh       initalize_ibm.sh      post_update.sh           resample_to_calc.sh  update_ibm.sh \
  inverse_functions.sh  post_update_single.sh    run_gromacs.sh       update_ibm_single.sh \
  do_shift_nb.pl            inverse.sh            postupd_pressure.sh      source_wrapper.pl    update_imc.sh \
  EX_to_CG.pl               linsolve.m            potential_to_gromacs.sh  source_wrapper.sh postupd_smooth.sh \
  solve_matlab.sh solve_octave.sh linsolve.octave merge_tables.pl smooth.pl postupd_splinesmooth.sh

inverse_DATA = MD5SUM csg_table

EXTRA_DIST = ${inverse_SCRIPTS} $(inverse_DATA)

