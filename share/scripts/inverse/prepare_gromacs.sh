#! /bin/bash

if [ "$1" = "--help" ]; then
  echo This script implemtents the function initialize
  echo for the Inverse Boltzmann Method
  echo Usage: ${0##*/} last_sim_dir
  echo USES: die cp run_or_exit grompp
  echo NEEDS:
  exit 0
fi

check_deps "$0"

[[ -z "$1" ]] && die "Missing argument for ${0##*/}"

cp ${1}/confout.gro ./conf.gro || die "${0##*/} cp ${1}/confout.gro ./conf.gro failed" 

run_or_exit grompp -n index.ndx 
