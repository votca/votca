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

# temporary compatibility issue

sed -ie 's/NaN/0.0/' $1.dpot.new
sed -ie 's/Inf/0.0/' $1.dpot.new 
