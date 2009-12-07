#! /bin/bash

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This script implemtents the function update
for the Inverse Boltzmann Method

Usage: ${0##*/} step_nr

USES:  die msg csg_get_property for_all do_external

NEEDS: cg.inverse.program 
EOF
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

msg "Calc rdf"
sim_prog="$(csg_get_property cg.inverse.program)" 
for_all non-bonded do_external rdf $sim_prog
for_all non-bonded do_external update ibm_single "${1}"
