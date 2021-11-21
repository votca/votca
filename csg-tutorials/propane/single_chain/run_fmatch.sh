#!/bin/bash

if [ ! -f confout.gro ]; then
  echo "confout.gro not found, please start run.sh first to run the simulation"
  exit 1
fi

csg_fmatch --top topol.tpr --trj traj.trr --cg propane.xml --options settings.xml

csg_call table integrate angle.force angle.pot.fm
csg_call table integrate bond.force bond.pot.fm

