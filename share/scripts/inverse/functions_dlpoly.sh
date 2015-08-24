#!/bin/bash
#
# Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
Useful functions for the generic simulation program:
EOF
sed -n 's/^\(.*\)([)] {[^#]*#\(.*\)$/* \1   -- \2/p' ${0}

  exit 0
fi

simulation_finish() { #checks if simulation is finished
  local nneeded npassed
  [[ ! -f "HISTORY" ]] && return 1
  [[ ! -s "OUTPUT"  ]] && return 1
  nneeded=$(awk '/selected number of timesteps/{print $5}' OUTPUT)
  npassed=$(awk '/run terminated after/{print $4}' OUTPUT)
  [[ $npassed -lt $nneeded ]] && return 1
  echo "DL_POLY simulation completed"
  traj=$(csg_get_property cg.inverse.dlpoly.traj)
  topol=$(csg_get_property cg.inverse.dlpoly.topol)
  critical touch $traj
  critical touch $topol
  return 0
}
export -f simulation_finish

checkpoint_exist() { #check if a checkpoint exists (REVIVE _and_ REVCON - both are needed!)
  #support for checkpoint
  local checkpoint check
  checkpoint="$(csg_get_property cg.inverse.dlpoly.checkpoint)"
  for check in $checkpoint; do
    [[ ! -f ${check} ]] && return 1
  done
  echo "DL_POLY checkpoint present (${checkpoint} found)"
  return 0
}
export -f checkpoint_exist

get_simulation_setting() { #gets parameter a parameter from the settings file (1st argument) from simulation setting file (not implemented)
  local sim_prog="$(csg_get_property cg.inverse.program)"
  die "${FUNCNAME[0]}: Not implemented for ${sim_prog} yet"
  return 1
}
export -f get_simulation_setting
