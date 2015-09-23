#! /bin/bash -e

grompp -v

mdrun -v

echo Calculating distributions
csg_stat --top topol.tpr --trj traj.trr --cg hexane.xml --nt 3 --options fmatch.xml

echo "Mapping confout.gro to get configuration for coarse-grained run"
csg_map --top topol.tpr --trj confout.gro --cg hexane.xml --out conf_cg.gro 
