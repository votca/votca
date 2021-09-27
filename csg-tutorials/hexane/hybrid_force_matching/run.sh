#! /bin/bash -e

if [[ ! -f ../atomistic/topol.tpr || ! -f ../atomistic/traj.trr ]]; then
  echo "Run atomistic simulation in ../atomistic first"
  exit 1
fi

echo "Rerun md trajectory with excluded bonded interactions"
gmx grompp -c ../atomistic/conf.gro -f ../atomistic/grompp.mdp
gmx mdrun -v -rerun ../atomistic/traj.trr


#equilibration time in Gromacs units (ps)
equi=200
echo equi = $equi

echo "Running force matching"
csg_fmatch --top topol.tpr --trj traj.trr --begin $equi  --options fmatch.xml --cg hexane.xml

#integrate force tables to get potential
csg_call table integrate A-A.force A-A.pot
csg_call table linearop A-A.pot A-A.pot -1 0

csg_call table integrate B-B.force B-B.pot
csg_call table linearop B-B.pot B-B.pot -1 0

csg_call table integrate A-B.force A-B.pot
csg_call table linearop A-B.pot A-B.pot -1 0

cp A-A.pot input_A-A.pot
cp A-B.pot input_A-B.pot
cp B-B.pot input_B-B.pot

csg_call --ia-type non-bonded --ia-name A-A --options fmatch.xml convert_potential gromacs --clean input_A-A.pot table_A_A.xvg
csg_call --ia-type non-bonded --ia-name A-B --options fmatch.xml convert_potential gromacs --clean input_A-B.pot table_A_B.xvg
csg_call --ia-type non-bonded --ia-name B-B --options fmatch.xml convert_potential gromacs --clean input_B-B.pot table_B_B.xvg

