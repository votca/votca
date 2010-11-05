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


dt=$(get_from_mdp dt "grompp.mdp")
equi_time="$(csg_get_property cg.inverse.gromacs.equi_time 0)"
first_frame="$(csg_get_property cg.inverse.gromacs.first_frame 0)"
nsteps=$(get_from_mdp nsteps "grompp.mdp")

name=$(csg_get_interaction_property name)


mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
[ -f "$mdp" ] || die "${0##*/}: gromacs mdp file '$mdp' not found"
adress_type=$(get_from_mdp adress_type "$mdp")

echo "Calculating $adress_type density"

#TODO implement property to read -sl from xml file
begin="$(awk -v dt=$dt -v frames=$first_frame -v eqtime=$equi_time 'BEGIN{print (eqtime > dt*frames ? eqtime : dt*frames) }')"


if [ $adress_type = "sphere" ]
then
dens_prog="g_rdf -n index.ndx -bin 0.01"
index_sel="DUM \n CG\n"
else
dens_prog="g_density -n index.ndx -sl 100 -d x"
index_sel=$name
fi

log "Running $dens_prog"
echo -e $index_sel | run_or_exit $dens_prog -b ${begin} -o dens.$name.xvg

	



