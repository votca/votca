#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script calcs the x-density for gromacs
   echo for the AdResS therm force
   echo Usage: ${0##*/}
   echo USES: get_from_mdp csg_get_property awk log run_or_exit g_energy csg_taillog die
   echo NEEDS: cg.inverse.gromacs.equi_time cg.inverse.gromacs.first_frame
   exit 0
fi


check_deps "$0"

nsteps=$(get_from_mdp nsteps)
dt=$(get_from_mdp dt)
equi_time="$(csg_get_property cg.inverse.gromacs.equi_time 0)"
first_frame="$(csg_get_property cg.inverse.gromacs.first_frame 0)"
name=$(csg_get_interaction_property name)

xstart="$(csg_get_property cg.tf.xstart)"
xstop="$(csg_get_property cg.tf.xstop)"
step="$(csg_get_interaction_property step)"
prefactor="$(csg_get_property cg.tf.prefactor)"
splinestep="$(csg_get_property cg.tf.splinesmoothstep)"

adressw="$(csg_get_property cg.tf.adressw)"
adressh="$(csg_get_property cg.tf.adressh)"
adressc="$(csg_get_property cg.tf.adressc)"


infile="dens_$name.xvg"
outfile="dens_$name\_smooth.xvg"

#rho_0="$(cat dens_$name.xvg | awk -f $CSGSCRIPTDIR/calc_non_hybrid_dens.awk -v adressc=$adressc adressw=$adressw adressh=$adressh)"
#log "Density in non hybrid zone $rho_0"

log "Symmetrizing density profile"
infile="dens_$name.xvg"
outfile="dens_$name\_symm.xvg"
run_or_exit do_external density symmetrize --infile $infile --outfile $infile --adressc $adressc

log "Smoothing density profile"
infile="dens_$name\_symm.xvg"
outfile="dens_$name\_smooth.xvg"
run_or_exit csg_resample --in $infile --out $outfile --grid $xstart:$step:$xstop --derivative $forcefile --spfit $xstart:$splinestep:$xstop

forcefile="thermforce_$name.xvg"
forcefile_smooth="thermforce_smooth_$name.xvg"

run_or_exit do_external table smooth_borders --infile $forcefile --outfile $forcefile_smooth --xstart $xstart --xstop $xstop  

run_or_exit do_external table integrate $forcefile_smooth ${name}.dpot.new $prefactor


