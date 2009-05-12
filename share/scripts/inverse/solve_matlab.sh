  #! /bin/bash

if [ "$1" = "--help" ]; then
   echo "This solves linear equation system from imc using matlab"
   echo "Usage: ${0##*/} <group>"
   echo "Needs:  matlab"
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing arguments"

r_min=$(sed -ne '/i[[:space:]]*$/p' CG-CG.pot.cur | sed -ne '1p' | awk '{print $1}')

sed -e "s/\$name/$1/" -e "s/\$r_min/$r_min/" $CSGINVERSE/linsolve.m > solve_$1.m

# kill the flags
awk '{print $1,$2}' $1.imc > ${1}_noflags.imc
matlab -r solve_$1 -nosplash -nodesktop
rm -f solve_$1.m

# add stupid flags
sed -ie '/^[#@]/d' $1.dpot.matlab
sed -ie '/[0-9]/s/$/ i/' $1.dpot.matlab

# copy flags
merge_tables="$($SOURCE_WRAPPER tools merge_tables)" || die "${0##*/}: $SOURCE_WRAPPER tools merge_tables failed"
run_or_exit $merge_tables --novalues $1.pot.cur $1.dpot.matlab $1.dpot.new

# temporary compatibility issue

sed -ie 's/NaN/0.0/' $1.dpot.new
sed -ie 's/Inf/0.0/' $1.dpot.new 
