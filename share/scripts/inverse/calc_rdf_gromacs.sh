#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script calcs the rdf for gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo Needs: g_rdf
   exit 0
fi

for exe in g_rdf; do
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
#echo equi = $equi

echo Running g_rdf for ${type1}-${type2}
echo -e "${type1}\n${type2}" | g_rdf -b ${equi} -n index.ndx -bin 0.01 -o rdf_${type1}_${type2}.xvg &> log_g_rdf_${type1}_${type2}
[[ $? -ne 0 ]] && echo Error at g_rdf ${type1}-${type2} > /dev/stderr && exit 1

