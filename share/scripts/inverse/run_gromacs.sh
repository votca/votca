#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script runs gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo USES: run_or_exit mdrun
   exit 0
fi

check_deps "$0"

run_or_exit mdrun
