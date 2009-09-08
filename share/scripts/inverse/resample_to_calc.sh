#! /bin/bash
if [ "$1" = "--help" ]; then
   echo This script resamples target distribution to grid spacing
   echo for calculations
   echo Usage: ${0##*/} target_directory
   echo USES:  die \$csg_get run_or_exit csg_resample
   echo NEEDS: min max step target name
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

min=$($csg_get min )
max=$($csg_get max )
step=$($csg_get step )
target=$($csg_get target)
name=$($csg_get name)

run_or_exit csg_resample --in ${target} --out ${1}/${name}.dist.tgt --grid ${min}:${step}:${max}
