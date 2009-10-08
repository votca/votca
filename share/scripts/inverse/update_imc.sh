#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function update
   echo for the Inverse Monte Carlo Method
   echo Usage: ${0##*/} step_nr
   echo USES: die csg_get_property msg is_done run_or_exit csg_stat \$CSGXMLFILE mark_done do_external sort 
   echo NEEDS: cg.cgmap cg.inverse.gromacs.topol cg.inverse.gromacs.traj cg.inverse.imc.solver cg.*.imc.group  cg.inverse.program
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

solver=$(csg_get_property cg.inverse.imc.solver)
sim_prog="$(csg_get_property cg.inverse.program)" 
do_external imc_stat $sim_prog

list_groups=$(csg_get_property 'cg.*.imc.group' | sort -u)
for group in "$list_groups"; do
  # currently this is a hack! need to create combined array
  msg "solving linear equations for $group"
  run_or_exit csg_imcrepack --in ${group} --out ${group}.packed
  do_external imcsolver $solver ${group}.packed ${group}.packed.sol
  run_or_exit csg_imcrepack --in ${group}.packed --unpack ${group}.packed.sol
done 

for_all non-bonded do_external imc purify 
