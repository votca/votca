#! /bin/bash

echo "Start simulaion"

grompp -v

mdrun -v

echo "Finished simulation"

echo "Calculating target RDFs"

echo "Calculating target RDF between urea-urea"
echo "1 1"|g_rdf -f traj.xtc -o UR-UR.dist.tgt -b 10000 -rdf mol_com -n

echo "Calculating target RDF between urea-water"
echo "1 2"|g_rdf -f traj.xtc -o UR-SOL.dist.tgt -b 10000 -rdf mol_com -n

echo "Calculating target RDF between water-water"
echo "2 2"|g_rdf -f traj.xtc -o SOL-SOL.dist.tgt -b 10000 -rdf mol_com -n -dt 10

echo "Mapping confout.gro to get configuration for coarse-grained run"
csg_map --top topol.tpr --trj confout.gro --cg "urea.xml;water.xml" --out conf_cg.gro

# end of run file

