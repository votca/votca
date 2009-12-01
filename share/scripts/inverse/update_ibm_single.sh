#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function update
   echo for the Inverse Boltzmann Method for a single pair
   echo Usage: ${0##*/} step_nr
   echo USES:  do_external die csg_get_interaction_property log run_or_exit awk
   echo NEEDS: name step min max inverse.do_potential
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

scheme=( $(csg_get_interaction_property inverse.do_potential 1) )
scheme_nr=$(( ( $1 - 1 ) % ${#scheme[@]} ))
name=$(csg_get_interaction_property name)

if [ "${scheme[$scheme_nr]}" = 1 ]; then
   log "Update potential ${name} : yes"
   #update ibm
   run_or_exit do_external update ibm_pot ${name}.dist.tgt ${name}.dist.new ${name}.pot.cur ${name}.dpot.tmp
   run_or_exit do_external dpot shift_nb ${name}.dpot.tmp ${name}.dpot.new
else
   log "Update potential ${name} : no"
   awk -v step=$(csg_get_interaction_property step) -v start=$(csg_get_interaction_property min) -v end=$(csg_get_interaction_property max) \
     'BEGIN{x=start;while(x<end+step){print x,0.0,"i";x+=step;}}' > ${name}.dpot.new \
      || die "${0##*/}: awk failed"
fi
