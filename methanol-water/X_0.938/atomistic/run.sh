#! /bin/bash -e

gmx grompp -v

gmx mdrun -v

equi=2000
echo equi = $equi

#determine number of threads nt to run csg_stat in parallel
nt="$(grep -c processor /proc/cpuinfo 2>/dev/null)" || nt=0
((nt++))

echo Calculating distributions
csg_stat --top topol.tpr --trj traj.trr --cg "methanol.xml;water.xml" --options settings.xml --nt $nt --begin $equi

echo "Mapping confout.gro to get configuration for coarse-grained run"
csg_map --top topol.tpr --trj confout.gro --cg "methanol.xml;water.xml" --out conf_cg.gro 
