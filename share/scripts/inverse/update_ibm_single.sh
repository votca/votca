#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function update
   echo for the Inverse Boltzmann Method for a single pair
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$source_wrapper, update_POT.pl
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

update_POT="$($SOURCE_WRAPPER --direct update_POT.pl)" || die "${0##*/}: $SOURCE_WRAPPER --direct update_POT.pl failed" 

scheme=( $($csg_get  .do_potential ) )
scheme_nr=$(( ( $1 - 1 ) % ${#scheme[@]} ))
type1=$($csg_get type1)
type2=$($csg_get type2)

if [ "${scheme[$scheme_nr]}" = 1 ]; then
   log "Update potential ${type1}-${type2} : yes"
   #update ibm
   run_or_exit ${update_POT} ${type1}_${type2}.dist.tgt ${type1}_${type2}.dist.new ${type1}_${type2}.dpot.new
else
   log "Update potential ${type1}-${type2} : no"
fi
