#! /bin/bash -e

echo "Start simulaion"

gmx grompp -v

gmx mdrun -v

echo "Finished simulation"

#equilibration time in Gromacs units (ps)
equi=10000
echo equi = $equi


#determine number of threads nt to run csg_stat in parallel
nt="$(grep -c processor /proc/cpuinfo 2>/dev/null)" || nt=0
((nt++))

echo "Calculating target RDFs"
csg_stat --nt $nt --top topol.tpr --cg "urea.xml;water.xml" --options non-bonded.xml --trj traj.xtc --begin $equi

echo "Mapping confout.gro to get configuration for coarse-grained run"
csg_map --top topol.tpr --trj confout.gro --cg "urea.xml;water.xml" --out conf_cg.gro

# end of run file

