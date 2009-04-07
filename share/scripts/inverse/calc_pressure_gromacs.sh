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
#20 % is warmup
equi=$(awk "BEGIN{print 0.2*$nsteps*$dt}")

log "Running g_energy"
echo "Pressure" | g_energy -b $equi >> $CSGLOG 2>&1 || die "${0##*/}: Error at running g_energy"
p_now=$(tail -30 $CSGLOG | awk '/^Pressure/{print $3}' ) || die "${0##*/}: awk failed"
echo ${p_now}
