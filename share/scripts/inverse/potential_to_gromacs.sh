#!/bin/bash

if [ "$1" = "--help" ]; then
  echo This is a wrapper to convert potential to gromacs
  echo Usage: ${0##*/} 
  echo Needs: SOURCE_WRAPPER, run_or_exit, wc 
  exit 0
fi

table_to_xvg="$($SOURCE_WRAPPER --direct table_to_xvg.pl)" || die "$SOURCE_WRAPPER --direct table_to_xvg.pl failed" 

name=$($csg_get name)
input="${name}.pot.cur" 
#gromacs want '_' !
output="table_$($csg_get type1)_$($csg_get type2).xvg" 
log "Convert $input to $output"

min=$($csg_get min)
max=$($csg_get max)
gromacs_bins="$(get_sim_property gromacs.table_bins)"

log echo running "csg_resample --in $input --out smooth_${input} --grid $min:$gromacs_bins:$max"

csg_resample --in $input --out smooth_${input} --grid $min:$gromacs_bins:$max || die "${0##*/}: csg_resample failed"


run_or_exit $table_to_xvg smooth_${input} $output 

