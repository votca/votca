#! /bin/bash

grompp -v

simnane=nvt
mdrun -append -cpi state.cpt &> log_mdrun 

echo Running g_energy
equi=2000
echo equi = $equi
echo -e "Pressure" | g_energy -b $equi &> log_g_energy
[[ $? -ne 0 ]] && echo Error at running g_energy && exit 1

echo Running g_rdf
echo -e "a O\nq" | make_ndx -f conf.gro &> log_make_ndx
echo -e "MeO\nMeO" | g_rdf -f topol.trr -b $equi -bin 0.01 -rdf mol_com -n index.ndx -o rdf_methanol.xvg &> log_g_rdf_CG_CG

echo "Mapping confout.gro to get configuration for coarse-grained run"
csg_map --top topol.tpr --trj confout.gro --cg methanol.xml --out conf_cg.gro 
