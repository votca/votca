
# BEGIN: these scripts have to be defined for each method
# method initialization
prepare tf prepare_generic.sh
prepare_single tf prepare_tf_single.sh
#init_single tf init_tf.sh
initstep tf initialize_step_generic.sh

# the update step
update tf update_tf.sh
update tf_single update_tf_single.sh

# add up the potential
# TODO: think about name add_pot, is other name maybe better
add_pot tf add_pot_generic.sh
# END: these scripts have to be defined for each method

calc thermforce calc_thermforce.sh

convert_potential xvg table_to_xvg.pl
tf apply_prefactor apply_prefactor.pl
# todo : can aplly_pref be the post update?
post_update tf dummy.sh


density gromacs calc_density_gromacs.sh
density symmetrize density_symmetrize.py 

table getsubset table_getsubset.py
table smooth_borders table_smooth_borders.py
