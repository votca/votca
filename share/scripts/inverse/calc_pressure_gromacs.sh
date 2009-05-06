#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script calcs the pressure for gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo Needs: g_energy
   exit 0
fi

for exe in g_energy awk; do
   if [ -z $(type -p $exe) ]; then
      echo Could not find $exe > /dev/stderr
      exit 1
   fi
done

if [ $(type -t get_from_mdp) != "function" ]; then
   echo Could not find function get_from_mdp > /dev/stderr
   exit 1
fi

nsteps=$(get_from_mdp nsteps)
dt=$(get_from_mdp dt)
equi_time="$(get_sim_property gromacs.equi_time)"
equi_time=${equi_time%\%}
equi=$(awk "BEGIN{print $equi_time/100*$nsteps*$dt}")

log "Running g_energy"
run_or_exit "echo Pressure | g_energy -b $equi"
p_now=$(tail -30 $CSGLOG | awk '/^Pressure/{print $3}' ) || die "${0##*/}: awk failed"
echo ${p_now}
