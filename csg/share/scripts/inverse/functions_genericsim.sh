#!/usr/bin/env bash
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
  local traj sim_prog
  sim_prog="$(csg_get_property cg.inverse.program)"
  traj="$(csg_get_property cg.inverse.$sim_prog.traj)"
  [[ -f $traj ]] && return 0 
  return 1
}
export -f simulation_finish

checkpoint_exist() { #check if a checkpoint exists (not implemented)
  #no support for checkpoints, yet !
  return 1
}
export -f checkpoint_exist

get_simulation_setting() { #gets parameter a parameter from the settings file (1st argument) from simulation setting file (not implemented)
  local sim_prog="$(csg_get_property cg.inverse.program)"
  die "${FUNCNAME[0]}: Not implemented for ${sim_prog} yet"
  return 1
}
export -f get_simulation_setting
