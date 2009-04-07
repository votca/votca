#!/bin/bash

if [ "$1" = "--help" ]; then
  echo This is a wrapper to convert potential to gromacs
  echo Usage: ${0##*/} 
  echo Needs: SOURCE_WRAPPER, run_or_exit, wc, spline 
  exit 0
fi

table_to_xvg="$($SOURCE_WRAPPER --direct table_to_xvg.pl)" || die "$SOURCE_WRAPPER --direct table_to_xvg.pl failed" 

for exe in spline wc; do
   if [ -z $(type -p $exe) ]; then
      die "${0##*/}: Could not find $exe"
   fi
done

type1=$($csg_get type1)
type2=$($csg_get type2)
input="${type1}_${type2}.pot.cur" 
output="table_${type1}_${type2}.xvg" 
log "Convert $input to $output"
n_lines=$(wc -l $input | awk '{print 5*($1-1)}')
log "Spline lines are $n_lines for ${type1}-${type2}"
spline -n $n_lines $input > smooth_${input} || die "${0##*/}: spline failed"
run_or_exit $table_to_xvg smooth_${input} $output 

