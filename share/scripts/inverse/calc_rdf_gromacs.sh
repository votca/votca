#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script calcs the rdf for gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo Needs: \$atoms, g_rdf
   exit 0
fi

for exe in make_ndx g_rdf; do
   if [ -z $(type -p $exe) ]; then
      echo Could not find $exe > /dev/stderr
      exit 1
   fi
done

if [ -z "$scriptdir" ]; then
   echo scriptdir not defined > /dev/stderr
   exit 1
fi

if [ -z "$last_dir" ]; then
   echo last_dir not defined > /dev/stderr
   exit 1
fi

if [ ${#atoms[@]} -eq 0 ]; then
   echo Could not get atomname for \$atoms > /dev/stderr
   exit 1
fi

for ((i=0;i<${#atoms[*]};i++)); do
   atom1=${atoms[$i]}
   for ((j=$i;j<${#atoms[*]};j++)); do
      atom2=${atoms[$j]}
      echo Running g_rdf for ${atom1}-${atom2}
      echo -e "${atom1}\n${atom2}" | g_rdf -b $equi -n index.ndx -bin 0.01 -o rdf_${atom1}_${atom2}.xvg &> log_g_rdf_${atom1}_${atom2}
      [[ $? -ne 0 ]] && echo Error at g_rdf ${atom1}-${atom2} && exit 1
   done
done