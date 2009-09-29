  #! /bin/bash

if [ "$1" = "--help" ]; then
   echo "This solves linear equation system from imc using matlab"
   echo "Usage: ${0##*/} <group> <outfile>"
   echo "USES:  die sed matlab rm \$SOURCE_WRAPPER run_or_exit \$CSGINVERSE"
   echo NEEDS:
   exit 0
fi

check_deps "$0"

[[ -n "$2" ]] || die "${0##*/}: Missing arguments"

# initialize & run the matlab file
sed -e "s/\$name_out/$2/" -e "s/\$name/$1/" $CSGINVERSE/linsolve.m > solve_$1.m
matlab -r solve_$1 -nosplash -nodesktop
rm -f solve_$1.m

# temporary compatibility issue
sed -ie 's/NaN/0.0/' $2
sed -ie 's/Inf/0.0/' $2 
