#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script runs gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo USES: run_or_exit mdrun
   exit 0
fi

check_deps "$0"

if use_mpi; then
  mpicmd=$(csg_get_property cg.inverse.mpi.cmd)
  run_or_exit $mpicmd mdrun -tablea table_tf.xvg
else
  run_or_exit mdrun -nt 8 -tablea table_tf.xvg
fi
