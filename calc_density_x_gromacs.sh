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


dt=$(get_from_mdp dt "grompp.mdp")
equi_time="$(csg_get_property cg.inverse.gromacs.equi_time 0)"
first_frame="$(csg_get_property cg.inverse.gromacs.first_frame 0)"
nsteps=$(get_from_mdp nsteps "grompp.mdp")

name=$(csg_get_interaction_property name)


mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
[ -f "$mdp" ] || die "${0##*/}: gromacs mdp file '$mdp' not found"
adress_type=$(get_from_mdp adress_type "$mdp")

echo "Adress type: $adress_type"


begin="$(awk -v dt=$dt -v frames=$first_frame -v eqtime=$equi_time 'BEGIN{print (eqtime > dt*frames ? eqtime : dt*frames) }')"


densigroup="$(csg_get_interaction_property inverse.gromacs.density_group "UNDEF")"

g_densopt="$(csg_get_property --allow-empty cg.inverse.gromacs.g_density_options "")"


if [ $densigroup = "UNDEF" ]
then
index_sel=$name
echo "Calculating density for $name"
else
index_sel=$densigroup
fi

if [ $adress_type = "sphere" ]
then
run_or_exit csg_spheredens --trj traj.trr --top topol.tpr --bin 0.01 --out dens.$name.xvg --begin ${begin}
#dens_prog="g_rdf -n index.ndx -bin 0.01"
#index_sel="DUM \n CG\n"
else
dens_prog="g_density -n index.ndx -d x $g_densopt"
echo "Running $dens_prog"
echo -e $index_sel | run_or_exit $dens_prog -b ${begin} -o dens.$name.xvg
fi




	



