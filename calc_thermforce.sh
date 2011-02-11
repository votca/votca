#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script calcs the x-density for gromacs
   echo for the AdResS therm force
   echo Usage: ${0##*/}
   echo USES: get_from_mdp csg_get_property awk run_or_exit g_energy csg_taillog die
   echo NEEDS: cg.inverse.gromacs.equi_time cg.inverse.gromacs.first_frame
   exit 0
fi


check_deps "$0"

name=$(csg_get_interaction_property name)


splinedelta="$(csg_get_property cg.tf.splinedelta)"
step="$(csg_get_interaction_property step)"
prefactor="$(csg_get_property cg.tf.prefactor)"


cg_prefactor="$(csg_get_property --allow-empty cg.tf.cg_prefactor)"


splinestep="$(csg_get_property cg.tf.splinesmoothstep)"

adressw="$(csg_get_property cg.tf.adressw)"
adressh="$(csg_get_property cg.tf.adressh)"
adressc="$(csg_get_property cg.tf.adressc)"

# ../ in the next line is dirty but in first step gromppp is not yet copied to step_000
mdp="../$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
[ -f "$mdp" ] || die "${0##*/}: gromacs mdp file '$mdp' not found"
adress_type=$(get_from_mdp adress_type "$mdp")

infile="dens.${name}.xvg"
outfile="dens.${name}.smooth.xvg"

xstart=$(echo "scale=8; $adressc+$adressw" | bc)
xstop=$(echo "scale=8; $adressc+$adressw+$adressh" | bc)

if [ ! $adress_type = "sphere" ]
then
infile="dens.${name}.xvg"
outfile="dens.${name}.symm.xvg"
run_or_exit do_external density symmetrize --infile $infile --outfile $outfile --adressc $adressc
infile="dens.${name}.symm.xvg"
#note : in the spehere case (no symmetrizing necessary) infile stays dens.${name}.xvg, so this gets used for next step
fi

outfile="dens.${name}.smooth.xvg"

forcefile="thermforce.${name}.xvg"
forcefile_pref="thermforce.wpref.${name}.xvg"
forcefile_smooth="thermforce.smooth.${name}.xvg"

spxstart=$(echo "scale=8; $adressc+$adressw-$splinedelta" | bc)
spxstop=$(echo "scale=8; $adressc+$adressw+$adressh+$splinedelta" | bc)

comment="$(get_table_comment)"

run_or_exit csg_resample --type cubic --in $infile --out $outfile --grid $spxstart:$step:$spxstop --derivative $forcefile --fitgrid $spxstart:$splinestep:$spxstop --comment "$comment"

if [ -z "$cg_prefactor" ];then
       echo "Using fixed prefactor $prefactor "	
       run_or_exit do_external tf apply_prefactor $forcefile $forcefile_pref $prefactor
else
       echo "Using linear interpolation of prefactors. Ex. pref: $prefactor CG. pref : $cg_prefactor"
       run_or_exit do_external tf apply_prefactor $forcefile $forcefile_pref $prefactor $cg_prefactor
fi

run_or_exit do_external table smooth_borders --infile $forcefile_pref --outfile $forcefile_smooth --xstart $xstart --xstop $xstop  

run_or_exit do_external table integrate $forcefile_smooth ${name}.dpot.new


#todo: make this less dirty:
touch $name.dist.new 


