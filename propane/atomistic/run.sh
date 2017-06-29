#! /bin/bash -e

gmx grompp -v

gmx mdrun -v

echo Running Extract_Energies.sh to extract all thermodynamic quantities from ener.edr
#equilibration time in Gromacs units (ps)
equi=200
echo equi = $equi
./Extract_Energies.sh "$equi" "Bond" "Angle" "Ryckaert-Bell."

#determine number of threads nt to run csg_stat in parallel
nt="$(grep -c processor /proc/cpuinfo 2>/dev/null)" || nt=0
((nt++))

echo Calculating distributions
csg_stat --top topol.tpr --trj traj.trr --cg propane.xml --nt $nt --options fmatch.xml --begin $equi

echo "Mapping confout.gro to get configuration for coarse-grained run"
csg_map --top topol.tpr --trj confout.gro --cg propane.xml --out conf_cg.gro 
