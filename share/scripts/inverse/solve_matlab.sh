  #! /bin/bash

if [ "$1" = "--help" ]; then
   echo "This solves linear equation system from imc using matlab"
   echo "Usage: ${0##*/} <group> <outfile>"
   echo USES:  die sed matlab rm run_or_exit \$CSGINVERSE mv
   echo NEEDS:
   exit 0
fi

check_deps "$0"

[[ -n "$2" ]] || die "${0##*/}: Missing arguments"

# initialize & run the matlab file
sed -e "s/\$name_out/$2/" -e "s/\$name/$1/" $CSGINVERSE/linsolve.m > solve_$1.m || die "${0##*/}: sed failed"

#matlab does not like -_. etc in filenames
run_or_exit mv solve_$1.m solve.m
run_or_exit matlab -r solve -nosplash -nodesktop
rm -f solve.m

# temporary compatibility issue
[[ -f "$2" ]] || die "Matlab failed"
run_or_exit sed -ie 's/NaN/0.0/' $2
run_or_exit sed -ie 's/Inf/0.0/' $2
