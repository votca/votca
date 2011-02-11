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

step_nr=$(get_current_step_nr)
scheme=( $(csg_get_interaction_property inverse.do_tf 1) )
scheme_nr=$(( ($step_nr - 1 ) % ${#scheme[@]} ))
name=$(csg_get_interaction_property name)

if [ "${scheme[$scheme_nr]}" = 1 ]; then
   echo "Update tf ${name} : yes"
#update ibm
    msg "Calc density"
    sim_prog="$(csg_get_property cg.inverse.program)" 
    do_external density $sim_prog
    do_external calc thermforce
else
   echo "Update tf ${name} : no"
   min=$(csg_get_interaction_property min)
   step=$(csg_get_interaction_property step)
   max=$(csg_get_interaction_property max)
   do_external table dummy ${min}:${step}:${max} ${name}.dpot.new
fi




