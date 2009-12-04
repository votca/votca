#! /bin/bash

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This script implemtents the function update
for the Inverse Monte Carlo Method

Usage: ${0##*/} step_nr

USES: die csg_get_property msg run_or_exit do_external sort for_all

NEEDS: cg.inverse.imc.solver cg.*.inverse.imc.group  cg.inverse.program
EOF
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

solver=$(csg_get_property cg.inverse.imc.solver)
sim_prog="$(csg_get_property cg.inverse.program)" 
do_external imc_stat $sim_prog

list_groups=$(csg_get_property 'cg.*.inverse.imc.group' | sort -u)
for group in "$list_groups"; do
  # currently this is a hack! need to create combined array
  msg "solving linear equations for $group"
  run_or_exit csg_imcrepack --in ${group} --out ${group}.packed
  do_external imcsolver $solver ${group}.packed ${group}.packed.sol
  run_or_exit csg_imcrepack --in ${group}.packed --unpack ${group}.packed.sol
done 

for_all non-bonded do_external imc purify 
