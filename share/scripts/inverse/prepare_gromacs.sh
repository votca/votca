#! /bin/bash

if [ "$1" = "--help" ]; then
  echo This script implemtents the function initialize
  echo for the Inverse Boltzmann Method
  echo Usage: ${0##*/} $last_sim_dir
  echo Needs: \$atoms, run_or_exit, make_ndx, grompp,sort,uniq
  exit 0
fi

if [ -z "$1" ]; then
  echo Missing argument for ${0##*/} > /dev/stderr
  exit 1
fi

for exe in make_ndx grompp; do
   if [ -z $(type -p $exe) ]; then
      echo Could not find $exe > /dev/stderr
      exit 1
   fi
done

cp ${1}/confout.gro ./conf.gro || die "${0##*/} cp ${1}/confout.gro ./conf.gro failed" 

#atoms="$(for_all non-bonded 'echo $($csg_get type1)'; for_all non-bonded 'echo #$($csg_get type2)' 2>&1 )" || die "for_all non-bonded  failed"
#atoms="$(echo "${atoms}" | sort | uniq )" || die "sort uniq failed"
#[[ -n "${atoms}" ]] || die "no atoms found"
#log "${0##*/}: Found atoms $atoms"
#run_or_exit "echo -e \"t ${atoms}\\nq\" | make_ndx -f conf.gro"
run_or_exit grompp -n index.ndx 
