#!/bin/bash -e

if [ ! -f confout.gro ]; then
  echo "confout.gro not found, please start run.sh first to run the simulation"
  exit 1
fi

csg_boltzmann --top topol.tpr --trj traj.trr --cg propane.xml <  boltzmann_cmds

