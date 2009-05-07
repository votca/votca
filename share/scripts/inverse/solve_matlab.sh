  #! /bin/bash

if [ "$1" = "--help" ]; then
   echo "This solves linear equation system from mc using matlab"
   echo "Usage: ${0##*/} <group>"
   echo "Needs:  matlab"
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing arguments"

sed -e "s/%NAME/$1/" $CSGSHARE/linsolve.m > solve_$1.m
# kill the flags
awk '{print $1,$2}' $1.imc > ${1}_noflags.imc
/sw/linux/suse/client/matlab/bin/matlab -arch=glnx86 -r solve_$1 -nosplash -nodesktop
rm -f solve_$1.m

paste $1.dpot.matlab $1.pot.cur | awk '{print $1,$2,$5}' > $1.dpot.new
# temporary compatibility issue

sed -ie 's/NaN/0.0/' $1.dpot.new
sed -ie 's/Inf/0.0/' $1.dpot.new 
