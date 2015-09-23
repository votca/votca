#!/bin/bash -e

MDRUN=mdrun
GROMPP=grompp
TRJCONV=trjconv

$GROMPP -c ../md/conf.gro
$MDRUN -v -rerun ../md/traj.trr
echo 0 | $TRJCONV -pbc mol -f traj.trr -o fixed_mols.trr -force

csg_fmatch --top topol.tpr --trj fixed_mols.trr --cg ../md/hexane.xml --options fmatch.xml

for i in A-A B-B A-B; do
  csg_call table integrate $i.force $i.pot
  csg_call table linearop $i.pot $i.pot -1 0 
done


