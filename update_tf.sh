#! /bin/bash
if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script implemtents the function update
for the Inverse Boltzmann Method

Usage: ${0##*/}

USES:  die msg csg_get_property for_all do_external

NEEDS: cg.inverse.program
EOF
   exit 0
fi

check_deps "$0"

check_deps "$0"



msg "Calc density"
sim_prog="$(csg_get_property cg.inverse.program)" 
for_all non-bonded do_external density $sim_prog
for_all non-bonded do_external calc thermforce

