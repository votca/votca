#! /bin/bash -e

grompp -v

mdrun -v

equi=2000
echo equi = $equi

echo Calculating distributions
csg_stat --top topol.tpr --trj traj.trr --cg "methanol.xml;water.xml" --options settings.xml --nt 3 --begin $equi

echo "Mapping confout.gro to get configuration for coarse-grained run"
csg_map --top topol.tpr --trj confout.gro --cg "methanol.xml;water.xml" --out conf_cg.gro 
