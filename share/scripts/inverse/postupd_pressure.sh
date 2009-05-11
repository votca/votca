#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the pressure update
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$SOURCE_WRAPPER
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

sim_prog=$(get_sim_property program)
p_now="$(do_external pressure $sim_prog)" 
p_target="$(get_sim_property p_target)"  
name=$($csg_get name)

log "New pressure $p_now"
log "Target pressure was $p_target"
  
pscheme=( $($csg_get .do_pressure ) )
pscheme_nr=$(( ( $1 - 1 ) % ${#pscheme[@]} ))

if [ "${pscheme[$pscheme_nr]}" = 1 ]; then
   log "Update pressure ${name} : yes"
   #calc pressure correction
   pressure_cor=$($SOURCE_WRAPPER pressure correction) || die "${0##*/}:$SOURCE_WRAPPER pressure correction failed"
   run_or_exit $pressure_cor $p_target $p_now pressure_cor.d 

   add_POT=$($SOURCE_WRAPPER add pot) || die "${0##*/}: $SOURCE_WRAPPER add pot failed" 

   run_or_exit $add_POT pressure_cor.d ${name}.dpot.cur ${name}.dpot.new
else
   log "Update pressure ${name} : no"
   run_or_exit cp ${name}.dpot.cur ${name}.dpot.new
fi
