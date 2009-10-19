#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script calcs the pressure for gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo USES: get_from_mdp csg_get_property awk log run_or_exit g_energy csg_taillog die
   echo NEEDS: cg.inverse.gromacs.equi_time
   exit 0
fi

check_deps "$0"

nsteps=$(get_from_mdp nsteps)
dt=$(get_from_mdp dt)
equi_time="$(csg_get_property cg.inverse.gromacs.equi_time 0)"
first_frame="$(csg_get_property cg.inverse.gromacs.first_frame 0)"

begin="$(awk -v dt=$dt -v frames=$first_frame -v eqtime=$equi_time 'BEGIN{print (eqtime > dt*frames ? eqtime : dt*frames) }')"

log "Running g_energy"
run_or_exit "echo Pressure | g_energy -b ${begin}"
p_now=$(csg_taillog -30 | awk '/^Pressure/{print $3}' ) || die "${0##*/}: awk failed"
echo ${p_now}
