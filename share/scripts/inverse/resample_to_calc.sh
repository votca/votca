#! /bin/bash
if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This script resamples target distribution to grid spacing
for calculations

Usage: ${0##*/} target_directory

USES:  die csg_get_interaction_property run_or_exit csg_resample

NEEDS: min max step inverse.target name
EOF
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

min=$(csg_get_interaction_property min )
max=$(csg_get_interaction_property max )
step=$(csg_get_interaction_property step )
target=$(csg_get_interaction_property inverse.target)
name=$(csg_get_interaction_property name)

run_or_exit csg_resample --in ${1}/${target} --out ${name}.dist.tgt --grid ${min}:${step}:${max}
