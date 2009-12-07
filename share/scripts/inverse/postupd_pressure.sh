#! /bin/bash

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This script implemtents the pressure update

Usage: ${0##*/} step_nr

USES:  die csg_get_property do_external csg_get_interaction_property log run_or_exit cp

NEEDS: cg.inverse.program name

OPTIONAL: inverse.post_update_options.pressure.type inverse.post_update_options.pressure.do
EOF
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

sim_prog="$(csg_get_property cg.inverse.program)"
name=$(csg_get_interaction_property name)

p_now="$(do_external pressure $sim_prog)" || die "${0##*/}: do_external pressure $sim_prog failed"
[ -z "$p_now" ] && die "${0##*/}: Could not get pressure from simulation"
log "New pressure $p_now"

ptype="$(csg_get_interaction_property inverse.post_update_options.pressure.type simple)"
pscheme=( $(csg_get_interaction_property inverse.post_update_options.pressure.do 1 ) )
pscheme_nr=$(( ( $1 - 1 ) % ${#pscheme[@]} ))

if [ "${pscheme[$pscheme_nr]}" = 1 ]; then
   log "Apply ${ptype} pressure correction for interaction ${name}"
   run_or_exit do_external pressure_cor $ptype $p_now pressure_cor.d 
   run_or_exit do_external table add pressure_cor.d ${name}.dpot.cur ${name}.dpot.new
else
   log "NO pressure correction for interaction ${name}"
   run_or_exit cp ${name}.dpot.cur ${name}.dpot.new
fi
