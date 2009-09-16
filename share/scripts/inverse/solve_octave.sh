  #! /bin/bash

if [ "$1" = "--help" ]; then
   echo "This solves linear equation system from imc using octave"
   echo "Usage: ${0##*/} <group> <copy-flags-from>"
   echo "if <copy-flags-from> is specified, flags are copied from that file"
   echo "otherwise every flag will be set to i"
   echo USES:  die sed awk octave \$SOURCE_WRAPPER run_or_exit rm \$CSGINVERSE
   echo NEEDS:
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing arguments"

r_min=$(sed -ne '/i[[:space:]]*$/p' CG-CG.pot.cur | sed -ne '1p' | awk '{print $1}')

sed -e "s/\$name/$1/" -e "s/\$r_min/$r_min/" $CSGINVERSE/linsolve.octave > solve_$1.octave

# kill the flags
awk '{print $1,$2}' $1.imc > ${1}_noflags.imc
octave solve_$1.octave
rm -f solve_$1.octave

# temporary compatibility issue
sed -ie 's/NaN/0.0/' $1.dpot.octave
sed -ie 's/Inf/0.0/' $1.dpot.octave

# add stupid flags
sed -ie '/^[#@]/d' $1.dpot.octave
sed -ie '/[0-9]/s/$/ i/' $1.dpot.octave

# copy flags
merge_tables="$($SOURCE_WRAPPER table merge)" || die "${0##*/}: $SOURCE_WRAPPER table merge failed"

if [ -n "$2" ]; then
    run_or_exit $merge_tables --novalues $2 $1.dpot.octave $1.dpot.new
fi



