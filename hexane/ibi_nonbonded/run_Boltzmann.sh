#!/bin/bash -e

#calculate bonded potentials with csg_boltzmann

csg_boltzmann --top ../atomistic/topol.tpr --trj ../atomistic/traj.trr --cg ../atomistic/hexane.xml < boltzmann_cmds

#smooth bonded potentials
csg_call --sloppy-tables table smooth bond.pot.ib input_bond.pot
csg_call --sloppy-tables table smooth angle.pot.ib input_angle.pot

#convert bonded potentials to GROMACS tables
if [ -d table ]; then
  rm -r table
fi
mkdir table

csg_call --ia-type bond --ia-name bond --options bond.xml convert_potential gromacs --clean input_bond.pot ./table/table_b1.xvg
csg_call --ia-type angle --ia-name angle --options angle.xml convert_potential gromacs --clean input_angle.pot ./table/table_a1.xvg
rm input_angle.pot
