#! /bin/bash

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This script implemtents the function initialize
for the Inverse Boltzmann Method

Usage: ${0##*/} last_sim_dir

USES: die cp run_or_exit grompp

NEEDS:
EOF
  exit 0
fi

check_deps "$0"

[[ -z "$1" ]] && die "Missing argument for ${0##*/}"

cp ${1}/confout.gro ./conf.gro || die "${0##*/} cp ${1}/confout.gro ./conf.gro failed" 

run_or_exit grompp -n index.ndx 
