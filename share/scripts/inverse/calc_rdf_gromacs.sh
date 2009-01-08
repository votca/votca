#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script calcs the rdf for gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo Needs: \$atoms, g_rdf
   exit 0
fi

for exe in g_rdf; do
   if [ -z $(type -p $exe) ]; then
      echo Could not find $exe > /dev/stderr
      exit 1
   fi
done

if [ ${#atoms[@]} -eq 0 ]; then
   echo Could not get atomname for \$atoms > /dev/stderr
   exit 1
fi

get_from_mdp() {
   [[ -n "$1" ]] || { echo What?; exit 1;}
   sed -n -e "s#[[:space:]]*$1[[:space:]]*=[[:space:]]*\(.*\)\$#\1#p" grompp.mdp | sed -e 's#;.*##'
}

nsteps=$(get_from_mdp nsteps)
dt=$(get_from_mdp dt)
#20 % is warmup
equi=$(awk "BEGIN{print 0.2*$nsteps*$dt}")
echo equi = $equi

for ((a1=0;a1<${#atoms[*]};a1++)); do
   atom1=${atoms[$a1]}
   for ((a2=$a1;a2<${#atoms[*]};a2++)); do
      atom2=${atoms[$a2]}
      echo Running g_rdf for ${atom1}-${atom2}
      echo -e "${atom1}\n${atom2}" | g_rdf -b $equi -n index.ndx -bin 0.01 -o rdf_${atom1}_${atom2}.xvg &> log_g_rdf_${atom1}_${atom2}
      [[ $? -ne 0 ]] && echo Error at g_rdf ${atom1}-${atom2} && exit 1
   done
done

