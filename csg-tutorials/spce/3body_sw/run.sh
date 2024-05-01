#! /usr/bin/env bash -e

#equilibration time in Gromacs units (ps)
equi=200
echo equi = $equi

if [[ ! -f ../atomistic/topol.tpr || ! -f ../atomistic/traj.trr ]]; then
  echo "Run atomistic simulation in ../atomistic first"
  exit 1
fi


#echo "Running force matching"
csg_fmatch --top ../atomistic/topol.tpr --trj ../atomistic/traj.trr --begin $equi  --options fmatch.xml --cg water.xml

#integrate force table to get potential
csg_call table switch_border CG-CG.force CG-CG.force_switched 1.0
csg_call table integrate CG-CG.force_switched CG-CG.pot
csg_call table linearop CG-CG.pot CG-CG.pot -1 0

#copy CG-CG.pot and CG-CG-CG.pot to new files to prevent deletion by clean command
cp CG-CG.pot input.pot
cp CG-CG-CG.pot input_angular.pot

#convert to lammps potentials
csg_call --options convert_tables.xml --ia-name CG-CG --ia-type non-bonded convert_potential lammps --clean input.pot table_CG_CG.txt

csg_call --options convert_tables.xml --ia-name CG-CG-CG --ia-type angle convert_potential lammps --clean --no-shift input_angular.pot table_CG_CG_CG.txt

#run the LAMMPS simulation (needs a current LAMMPS version compiled with the user pair_style sw/table)
lmp < spce.in > spce.out

#calculate the RDF and angular distribution function
csg_stat --options calculate_distributions.xml --top traj.dump --trj traj.dump
