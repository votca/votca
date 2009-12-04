#! /bin/bash

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This script runs gromacs
for the Inverse Boltzmann Method

Usage: ${0##*/}

USES: run_or_exit mdrun
EOF
   exit 0
fi

check_deps "$0"

run_or_exit mdrun
