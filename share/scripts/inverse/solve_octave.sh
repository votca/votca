  #! /bin/bash

if [ "$1" = "--help" ]; then
   echo "This solves linear equation system from imc using matlab"
   echo "Usage: ${0##*/} <group> <outfile>"
   echo "USES:  die sed octave rm \$SOURCE_WRAPPER run_or_exit \$CSGINVERSE"
   echo NEEDS:
   exit 0
fi

check_deps "$0"

[[ -n "$2" ]] || die "${0##*/}: Missing arguments"

# initialize & run the matlab file
sed -e "s/\$name_out/$2/"  -e "s/\$name/$1/" $CSGINVERSE/linsolve.octave > solve_$1.octave || die "${0##*/}: sed failed"
run_or_exit octave solve_$1.octave
#rm -f solve_$1.octave

[[ -f "$2" ]] || die "Octave failed"
# temporary compatibility issue
run_or_exit sed -ie 's/NaN/0.0/' $2
run_or_exit sed -ie 's/Inf/0.0/' $2

