#! /bin/bash

if [ "$1" = "--help" ]; then
   echo "This script implemtents smoothing of the potential update (.dpot)"
   echo "Usage: ${0##*/} step_nr"
   echo "USES: die \$csg_get mktemp sed awk csg_resample"
   echo "NEEDS: name min max step "
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

name=$($csg_get name)
min=$($csg_get min)
max=$($csg_get max)
step=$($csg_get step)

tmpfile=$(mktemp ${name}.XXX) || die "mktemp failed"
  
sed -ne '/i[[:space:]]*$/p' CG-CG.dpot.cur > $tmpfile
spmin=$(sed -ne '1p' $tmpfile | awk '{print $1}')
spmax=$(sed -ne '$p' $tmpfile | awk '{print $1}')
spstep=$($csg_get post_update_options.splinesmooth.step)

csg_resample --in $tmpfile --out $name.dpot.new --grid $min:$step:$max --spfit $spmin:$spstep:$spmax

