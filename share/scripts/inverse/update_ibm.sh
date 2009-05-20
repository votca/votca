#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function update
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$SOURCE_WRAPPER
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

msg "Calc rdf"
sim_prog="$(csg_get_property inverse.program)" 
do_external rdf $sim_prog for_all non-bonded
do_external update ibm_single for_all non-bonded "${1}"

