#!/bin/bash

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This scripts cleans up the dpot tables for each interaction when using IMC

Usage: ${0##*/}

USES:  do_external run_or_exit csg_get_interaction_property csg_get_property log csg_resample

NEEDS: name min max step cg.inverse.kBT inverse.do_potential
EOF
   exit 0
fi

check_deps "$0"

name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)
step=$(csg_get_interaction_property step)
kBT=$(csg_get_property cg.inverse.kBT)
log "purifying dpot for $name"

run_or_exit csg_resample --in ${name}.dpot.imc --out ${name}.dpot.impure --grid ${min}:${step}:${max}

scheme=( $(csg_get_interaction_property inverse.do_potential 1) )
scheme_nr=$(( ( $1 - 1 ) % ${#scheme[@]} ))

if [ "${scheme[$scheme_nr]}" = 1 ]; then
  log "Update potential ${name} : yes"
  run_or_exit do_external table linearop --withflag o ${name}.dpot.impure ${name}.dpot.impure 0 0
  run_or_exit do_external table linearop --withflag i ${name}.dpot.impure ${name}.dpot.impure $kBT 0

  do_external dpot crop  ${name}.dpot.impure  ${name}.dpot.after_crop
  run_or_exit do_external dpot shift_nb ${name}.dpot.after_crop ${name}.dpot.new
else
  log "Update potential ${name} : no"
  run_or_exit do_external table linearop ${name}.dpot.impure ${name}.dpot.new 0 0
fi

