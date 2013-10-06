#!/bin/bash -e

MDRUN=mdrun
GROMPP=grompp
TRJCONV=trjconv

$GROMPP
$MDRUN -v
echo 0 | $TRJCONV -pbc mol -f traj.trr -o fixed_mols.trr -force

csg_fmatch --top topol.tpr --trj fixed_mols.trr --cg hexane.xml --options fmatch_full.xml

for i in A-A B-B A-B bond angle; do
  csg_call table integrate $i.force $i.pot
  csg_call table linearop $i.pot $i.pot -1 0 
done


