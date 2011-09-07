#!/bin/bash
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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
Useful functions for espresso:
EOF
sed -n 's/^\(.*\)([)] {[^#]*#\(.*\)$/* \1   -- \2/p' ${0}

echo
echo Used external packages: espresso
  exit 0
fi

#check performed when this file is sourced
esp_dir="$(csg_get_property --allow-empty cg.inverse.espresso.scriptdir)" || exit 1
if [[ -n ${esp_dir} ]]; then
  export ESPRESSO_SCRIPTS="${esp_dir}"
  [[ -d ${ESPRESSO_SCRIPTS} ]] || die "${BASH_SOURCE[0]}: cg.inverse.espresso.scriptdir ($ESPRESSO_SCRIPTS) is not a directory"
fi
unset esp_dir

simulation_finish() { #checks if simulation is finished
  local esp_success traj_esp espout
  esp_success="$(csg_get_property cg.inverse.espresso.success "success.esp")"
  traj_esp="$(csg_get_property cg.inverse.espresso.traj "top_traj.esp")"
  espout="$(csg_get_property cg.inverse.espresso.blockfile_out "confout.esp.gz")"
  [[ -f $esp_success && -f $traj_esp && -f $espout ]] && return 0 
  return 1
}
export -f simulation_finish

checkpoint_exist() { #check if a checkpoint exists
  #espresso has not support for checkpoints, yet !
  return 1
}
export -f checkpoint_exist

get_simulation_setting() { #gets parameter a parameter from the settings file (1st argument) from simulation setting file
  die "${FUNCNAME[0]}: Not implemented for Espresso yet"
  return 1
}
export -f get_simulation_setting
