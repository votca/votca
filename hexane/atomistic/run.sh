#! /bin/bash -e

gmx grompp -v

gmx mdrun -v

echo Running Extract_Energies.sh to extract all thermodynamic quantities from ener.edr
#equilibration time in Gromacs units (ps)
equi=200
echo equi = $equi
./Extract_Energies.sh $equi

echo Calculating distributions
csg_stat --top topol.tpr --trj traj.trr --cg hexane.xml --nt 3 --options settings.xml --begin $equi

echo "Mapping confout.gro to get the starting configurations for the coarse-grained runs of ibi"
csg_map --top topol.tpr --trj confout.gro --cg hexane.xml --out conf_cg.gro 
