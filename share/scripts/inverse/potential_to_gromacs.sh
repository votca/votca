#!/bin/bash

if [ "$1" = "--help" ]; then
  echo This is a wrapper to convert potential to gromacs
  echo Usage: ${0##*/} 
  echo USES: do_external csg_get_interaction_property log csg_get_property run_or_exit csg_resample
  echo NEEDS: name inverse.gromacs.table max cg.inverse.gromacs.table_bins
  exit 0
fi

check_deps "$0"

name=$(csg_get_interaction_property name)
input="${name}.pot.cur" 
#gromacs want '_' !
output="$(csg_get_interaction_property inverse.gromacs.table)" 
log "Convert $input to $output"

r_cut=$(csg_get_interaction_property max)
gromacs_bins="$(csg_get_property cg.inverse.gromacs.table_bins)"

run_or_exit csg_resample --in ${input} --out smooth_${input} --grid 0:${gromacs_bins}:${r_cut} 
run_or_exit do_external convert_potential xvg smooth_${input} ${output}
