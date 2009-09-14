#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function update
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/} step_nr
   echo USES:  die msg csg_get_property for_all do_external
   echo NEEDS: cg.inverse.program 
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

msg "Calc rdf"
sim_prog="$(csg_get_property cg.inverse.program)" 
for_all non-bonded do_external rdf $sim_prog
for_all non-bonded do_external update ibm_single "${1}"
