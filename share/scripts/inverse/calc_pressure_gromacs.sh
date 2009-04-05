#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script calcs the pressure for gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo Needs: g_energy
   exit 0
fi

for exe in g_energy; do
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

echo Running g_energy > /dev/stderr
echo "Pressure" | g_energy -b $equi &> log_g_energy
[[ $? -ne 0 ]] && echo Error at running g_energy > /dev/stderr && exit 1
p_now=$(awk '/^Pressure/{print $3}' log_g_energy)
echo ${p_now}
