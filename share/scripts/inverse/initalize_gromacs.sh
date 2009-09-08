#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This is a intialize stuff for gromacs
   echo Usage: ${0##*/}
   echo USES: cp die
   echo NEEDS:
   exit 0
fi

check_deps "$0"

cp ../conf.gro confout.gro || die "cp ../conf.gro confout.gro"
