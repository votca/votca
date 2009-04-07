#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script make all the post update with backup 
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$source_wrapper, add_POT.pl
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

sim_prog=$(get_sim_property program)
p_now="$(do_external pressure $sim_prog)" 
p_target="$(get_sim_property p_target)"  

log "New pressure $p_now"
log "Target pressure was $p_target"
  
#calc pressure correction
pressure_cor=$($SOURCE_WRAPPER --direct pressure_cor.pl) || die "${0##*/}:$SOURCE_WRAPPER --direct pressure_cor.pl failed"
run_or_exit $pressure_cor $p_target $p_now pressure_cor.d 

update_single=$($SOURCE_WRAPPER --direct post_update_single.sh) || die "${0##*/}: $SOURCE_WRAPPER --direct post_update_single.sh"
for_all non-bonded ${update_single} $1

