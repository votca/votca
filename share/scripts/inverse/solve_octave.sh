  #! /bin/bash

if [ "$1" = "--help" ]; then
   echo "This solves linear equation system from imc using octave"
   echo "Usage: ${0##*/} <group>"
   echo "Needs:  octave"
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing arguments"

r_min=$(sed -ne '/i[[:space:]]*$/p' CG-CG.pot.cur | sed -ne '1p' | awk '{print $1}')

sed -e "s/\$name/$1/" -e "s/\$r_min/$r_min/" $CSGSHARE/linsolve.octave > solve_$1.octave

# kill the flags
awk '{print $1,$2}' $1.imc > ${1}_noflags.imc
octave solve_$1.octave
rm -f solve_$1.octave

# add stupid flags
sed -ie '/^[#@]/d' $1.dpot.octave
sed -ie '/[0-9]/s/$/ i/' $1.dpot.octave

# copy flags
merge_tables="$($SOURCE_WRAPPER tools merge_tables)" || die "${0##*/}: $SOURCE_WRAPPER tools merge_tables failed"
$merge_tables --novalues $1.pot.cur $1.dpot.octave $1.dpot.new

# temporary compatibility issue

sed -ie 's/NaN/0.0/' $1.dpot.new
sed -ie 's/Inf/0.0/' $1.dpot.new 
