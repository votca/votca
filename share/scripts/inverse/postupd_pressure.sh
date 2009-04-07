#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the pressure update
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$source_wrapper, add_POT.pl
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

add_POT=$($SOURCE_WRAPPER --direct add_POT.pl) || die "${0##*/}: $SOURCE_WRAPPER --direct add_POT.pl failed" 

pscheme=( $($csg_get .do_pressure ) )
pscheme_nr=$(( ( $1 - 1 ) % ${#pscheme[@]} ))
type1=$($csg_get type1)
type2=$($csg_get type2)

if [ "${pscheme[$pscheme_nr]}" = 1 ]; then
   log "Update presuure ${type1}-${type2} : yes"
   run_or_exit $add_POT pressure_cor.d ${type1}_${type2}.dpot.cur ${type1}_${type2}.dpot.new
else
   log "Update pressure ${type1}-${type2} : no"
   run_or_exit cp ${type1}_${type2}.dpot.cur ${type1}_${type2}.dpot.new
fi
