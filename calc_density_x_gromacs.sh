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
nsteps=$(get_from_mdp nsteps "grommp.mdp")

name=$(csg_get_interaction_property name)


#TODO implement property to read -sl from xml file
begin="$(awk -v dt=$dt -v frames=$first_frame -v eqtime=$equi_time 'BEGIN{print (eqtime > dt*frames ? eqtime : dt*frames) }')"
log "Running g_density"
run_or_exit "echo ${name} | g_density -b ${begin} -d x -sl 500 -o dens_$name.xvg"



