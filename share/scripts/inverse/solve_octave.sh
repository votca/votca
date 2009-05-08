  #! /bin/bash

if [ "$1" = "--help" ]; then
   echo "This solves linear equation system from imc using octave"
   echo "Usage: ${0##*/} <group>"
   echo "Needs:  octave"
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing arguments"

sed -e "s/%NAME/$1/" $CSGSHARE/linsolve.octave > solve_$1.octave
# kill the flags
awk '{print $1,$2}' $1.imc > ${1}_noflags.imc
octave solve_$1.octave
rm -f solve_$1.octave

# copy flags
merge_tables="$($SOURCE_WRAPPER tools merge_tables)" || die "${0##*/}: $SOURCE_WRAPPER tools merge_tables failed"
run_or_exit $merge_tables --novalues $1.pot.cur $1.dpot.matlab $1.dpot.new

# temporary compatibility issue

sed -ie 's/NaN/0.0/' $1.dpot.new
sed -ie 's/Inf/0.0/' $1.dpot.new 
