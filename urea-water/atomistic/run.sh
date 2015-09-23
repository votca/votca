#! /bin/bash -e

echo "Start simulaion"

grompp -v

mdrun -v

echo "Finished simulation"

echo "Calculating target RDFs"
csg_stat --nt 5 --top topol.tpr --cg "urea.xml;water.xml" --options non-bonded.xml --trj traj.xtc --begin 10000

echo "Mapping confout.gro to get configuration for coarse-grained run"
csg_map --top topol.tpr --trj confout.gro --cg "urea.xml;water.xml" --out conf_cg.gro

# end of run file

