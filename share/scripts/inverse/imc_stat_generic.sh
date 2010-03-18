#! /bin/bash
#
# Copyright 2009 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script implemtents statistical analysis for the Inverse Monte Carlo Method
using generic csg tools

Usage: ${0##*/}

USES: msg run_or_exit mark_done csg_stat csg_get_property \$CSGXMLFILE is_done check_deps

NEEDS: cg.inverse.program cg.inverse.cgmap

OPTIONAL: cg.inverse.\$sim_prog.first_frame cg.inverse.\$sim_prog.equi_time cg.inverse.\$sim_prog.topol
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
cgmap=$(csg_get_property cg.inverse.cgmap)
topol=$(csg_get_property cg.inverse.$sim_prog.topol topol.tpr)
traj=$(csg_get_property cg.inverse.$sim_prog.traj traj.xtc)

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
