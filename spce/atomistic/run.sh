#! /bin/bash -e

gmx grompp -v

gmx mdrun -v

echo Running Extract_Energies.sh to extract all thermodynamic quantities from ener.edr
#equilibration time in Gromacs units (ps)
equi=200
echo equi = $equi
./Extract_Energies.sh $equi

#determine number of threads nt to run csg_stat in parallel
nt="$(grep -c processor /proc/cpuinfo 2>/dev/null)" || nt=0
((nt++))

echo Calculating distributions
csg_stat --top topol.tpr --trj traj.trr --cg water.xml --options settings.xml --nt $nt --begin $equi

echo "Mapping confout.gro to get the starting configuration for coarse-grained runs of ibi/imc/re"
csg_map --top topol.tpr --trj confout.gro --cg water.xml --out conf_cg.gro 
