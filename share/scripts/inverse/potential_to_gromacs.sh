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

name=$($csg_get name)
input="${name}.pot.cur" 
#gromacs want '_' !
output="table_$($csg_get type1)_$($csg_get type2).xvg" 
log "Convert $input to $output"

binsize=$($csg_get step)
gromacs_bins="$(get_sim_property gromacs.table_bins)"
factor=$(awk "BEGIN{print $binsize/$gromacs_bins}") || die "${0##*/}: awk failed"
log "Spline factor is $factor for ${name}"

n_lines=$(wc -l $input | awk -v factor=$factor '{print factor*($1-1)}')
log "Spline lines are $n_lines for ${name}"
spline -n $n_lines $input > smooth_${input} || die "${0##*/}: spline failed"
run_or_exit $table_to_xvg smooth_${input} $output 

