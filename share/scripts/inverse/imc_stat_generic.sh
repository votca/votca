#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents statistical analysis for the Inverse Monte Carlo Method
   echo using generic csg tools
   echo Usage: ${0##*/}
   echo USES: msg run_or_exit mark_done csg_stat
   echo NEEDS: cg.inverse.program cg.cgmap cg.inverse.sim_prog.topol cg.inverse.sim_prog.topol
   echo OPTIONAL: cg.inverse.$sim_prog.first_frame cg.inverse.$sim_prog.equi_time
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)" 
cgmap=$(csg_get_property cg.cgmap)
topol=$(csg_get_property cg.inverse.$sim_prog.topol)
traj=$(csg_get_property cg.inverse.$sim_prog.traj)

equi_time="$(csg_get_property cg.inverse.$sim_prog.equi_time 0)"
first_frame="$(csg_get_property cg.inverse.$sim_prog.first_frame 0)"

check_deps "$0"
msg "calculating statistics"
if is_done "imc_analysis"; then
  msg "IMC analysis is already done"
else
  msg "Running IMC analysis"
  run_or_exit csg_stat --do-imc --options $CSGXMLFILE --top $topol --trj $traj --cg $cgmap \
        --begin $equi_time --first-frame $first_frame
  mark_done "imc_analysis"
fi
