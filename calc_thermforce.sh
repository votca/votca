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

name=$(csg_get_interaction_property name)


splinedelta="$(csg_get_property cg.tf.splinedelta)"
step="$(csg_get_interaction_property step)"
prefactor="$(csg_get_property cg.tf.prefactor)"
splinestep="$(csg_get_property cg.tf.splinesmoothstep)"

adressw="$(csg_get_property cg.tf.adressw)"
adressh="$(csg_get_property cg.tf.adressh)"
adressc="$(csg_get_property cg.tf.adressc)"


infile="dens_$name.xvg"
outfile="dens_$name\_smooth.xvg"

xstart=$(echo "scale=8; $adressc+$adressw" | bc)
xstop=$(echo "scale=8; $adressc+$adressw+$adressh" | bc)

#rho_0="$(cat dens_$name.xvg | awk -f $CSGSCRIPTDIR/calc_non_hybrid_dens.awk -v adressc=$adressc adressw=$adressw adressh=$adressh)"
#log "Density in non hybrid zone $rho_0"


infile="dens_$name.xvg"
outfile="dens_$name\_symm.xvg"
run_or_exit do_external density symmetrize --infile $infile --outfile $outfile --adressc $adressc


infile="dens_$name\_symm.xvg"
outfile="dens_$name\_smooth.xvg"

forcefile="thermforce_$name.xvg"
forcefile_smooth="thermforce_smooth_$name.xvg"
spxstart=$(echo "scale=8; $adressc+$adressw-$splinedelta" | bc)
spxstop=$(echo "scale=8; $adressc+$adressw+$adressh+$splinedelta" | bc)
run_or_exit csg_resample --in $infile --out $outfile --grid $spxstart:$step:$spxstop --derivative $forcefile --spfit $spxstart:$splinestep:$spxstop

run_or_exit do_external table smooth_borders --infile $forcefile --outfile $forcefile_smooth --xstart $xstart --xstop $xstop  

run_or_exit do_external table integrate $forcefile_smooth ${name}.dpot.new $prefactor


