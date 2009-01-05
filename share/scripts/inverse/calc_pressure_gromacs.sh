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

get_from_mdp() {
   [[ -n "$1" ]] || { echo What?; exit 1;}
   sed -n -e "s#[[:space:]]*$1[[:space:]]*=[[:space:]]*\(.*\)\$#\1#p" grompp.mdp | sed -e 's#;.*##'
}

nsteps=$(get_from_mdp nsteps)
dt=$(get_from_mdp dt)
#20 % is warmup
equi=$(awk "BEGIN{print 0.2*$nsteps*$dt}")
echo equi = $equi

echo Running g_energy
echo "Pressure" | g_energy -b $equi &> log_g_energy
[[ $? -ne 0 ]] && echo Error at running g_energy && exit 1
p_now=$(awk '/^Pressure/{print $3}' log_g_energy)
