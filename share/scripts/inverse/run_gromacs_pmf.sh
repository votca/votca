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
This script runs gromacs
for the Inverse Boltzmann Method

Usage: ${0##*/}

USES: critical get_number_tasks csg_get_property check_deps

EOF
   exit 0
fi

tpr="$(csg_get_property cg.inverse.gromacs.topol "topol.tpr")"
[ -f "$tpr" ] || die "${0##*/}: gromacs tpr file '$tpr' not found"

mdrun="$(csg_get_property cg.inverse.gromacs.mdrun.bin "mdrun")"
[ -n "$(type -p $mdrun)" ] || die "${0##*/}: mdrun binary '$mdrun' not found"

confout="$(csg_get_property cg.inverse.gromacs.conf_out "confout.gro")"
opts="$(csg_get_property --allow-empty cg.inverse.gromacs.mdrun.opts)"

tasks=$(get_number_tasks)
mpicmd=$(csg_get_property --allow-empty cg.inverse.parallel.cmd)

# If in step_002, do mdrun -noappend
this_dir=$(get_current_step_dir --no-check)
pmf_step=$(basename $this_dir)
if [[ "$pmf_step" == "step_002" ]]; then
  # Rerun traj
  if [ -f "sim_done" ]; then
    rm sim_done
    traj=$(basename $(find . -maxdepth 1 -name traj.part*) )
    mpicmd="q2start -8tn rerun -f rerun_done --nompi"
    mdrun="mdrun -rerun ${traj} -noappend -cpi state.cpt -maxh 36 -e ener2.edr"
    opts="-pf pullf2.xvg"
  # Start/Continue run
  elif [ -f "rerun_done" ] || [ -f "prep_done" ]; then
    if [ -f "rerun_done" ]; then
      rm rerun_done
    elif [ -f "prep_done" ]; then
      rm prep_done
    fi
    traj=$(basename $(find . -maxdepth 1 -name traj.part*) )
    mpicmd="q2start -8tn run -f sim_done --nompi"
    mdrun="mdrun -noappend -cpi state.cpt -maxh 36 -e ener.edr"
    opts="-pf pullf.xvg"
  fi
fi

# Run gromacs
if [ $tasks -gt 1 ]; then
  critical $mpicmd $mdrun -s "${tpr}" -c "${confout}" ${opts}
else
  critical $mdrun -s "${tpr}" -c "${confout}" ${opts}
fi
