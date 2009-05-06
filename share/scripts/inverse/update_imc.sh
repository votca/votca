#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function update
   echo for the Inverse Monte Carlo Method
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$source_wrapper
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"


cgmap=$(get_property cgmap)
topol=$(get_sim_property gromacs.topol)
traj=$(get_sim_property gromacs.traj)
solver=$(get_sim_property imc.solver)
solver_bin=$($SOURCE_WRAPPER imcsolver $solver) || die "${0##*/}: $SOURCE_WRAPPER imcsolver $solver failed"


msg "calculating statistics"
run_or_exit csg_imc --options $CSGXMLFILE --top $topol --trj $traj --cg $cgmap

list_groups=$(csg_property --short --file $CSGXMLFILE --path "cg.*.imc.group" --print . | sort -u)
for group in "$list_groups"; do
  # currently this is a hack! need to create combined array
  msg "solving linear equations for $group"
  run_or_exit $solver_bin $group
done 
