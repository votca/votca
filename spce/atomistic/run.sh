#! /bin/bash -e

grompp -v

mdrun -v

echo Running g_energy
equi=2000
echo equi = $equi
echo -e "Pressure" | g_energy -b $equi

echo Calculating distributions
csg_stat --top topol.tpr --trj traj.trr --cg water.xml --options fmatch.xml --nt 3 --begin $equi

echo "Mapping confout.gro to get configuration for coarse-grained run"
csg_map --top topol.tpr --trj confout.gro --cg water.xml --out conf_cg.gro 
