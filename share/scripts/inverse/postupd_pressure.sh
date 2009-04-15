#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the pressure update
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$SOURCE_WRAPPER, add_POT.pl
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

add_POT=$($SOURCE_WRAPPER --direct add_POT.pl) || die "${0##*/}: $SOURCE_WRAPPER --direct add_POT.pl failed" 

pscheme=( $($csg_get .do_pressure ) )
pscheme_nr=$(( ( $1 - 1 ) % ${#pscheme[@]} ))
name=$($csg_get name)

if [ "${pscheme[$pscheme_nr]}" = 1 ]; then
   log "Update presuure ${name} : yes"
   run_or_exit $add_POT pressure_cor.d ${name}.dpot.cur ${name}.dpot.new
else
   log "Update pressure ${name} : no"
   run_or_exit cp ${name}.dpot.cur ${name}.dpot.new
fi
