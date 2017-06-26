#! /bin/bash -e

gmx grompp -v

gmx mdrun -v

echo Running Extract_Energies.sh to extract all thermodynamic quantities from ener.edr
#equilibration time in Gromacs units (ps)
equi=200
echo equi = $equi
./Extract_Energies.sh $equi


echo Calculating distributions
csg_stat --top topol.tpr --trj traj.trr --cg water.xml --options settings.xml --nt 3 --begin $equi

echo "Mapping confout.gro to get the starting configuration for coarse-grained runs of ibi/imc/re"
csg_map --top topol.tpr --trj confout.gro --cg water.xml --out conf_cg.gro 
