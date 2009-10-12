#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the pressure update
   echo Usage: ${0##*/} step_nr
   echo USES:  die csg_get_property do_external csg_get_interaction_property log run_or_exit cp
   echo NEEDS: cg.inverse.program cg.inverse.p_target name do_pressure
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

sim_prog=$(csg_get_property cg.inverse.program)
p_now="$(do_external pressure $sim_prog)" 
p_target="$(csg_get_property cg.inverse.p_target)"  
name=$(csg_get_interaction_property name)

log "New pressure $p_now"
log "Target pressure was $p_target"
  
pscheme=( $(csg_get_interaction_property do_pressure ) )
pscheme_nr=$(( ( $1 - 1 ) % ${#pscheme[@]} ))

if [ "${pscheme[$pscheme_nr]}" = 1 ]; then
   log "Update pressure ${name} : yes"
   #calc pressure correction
   run_or_exit do_external pressure correction $p_target $p_now pressure_cor.d 
   run_or_exit do_external table add pressure_cor.d ${name}.dpot.cur ${name}.dpot.new
else
   log "Update pressure ${name} : no"
   run_or_exit cp ${name}.dpot.cur ${name}.dpot.new
fi
