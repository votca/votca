#! /bin/bash
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
This script initializes an iteration for the generic simulation program 

Usage: ${0##*/}
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
from=$(csg_get_property cg.inverse.initial_configuration)
conf="$(csg_get_property --allow-empty cg.inverse.$sim_prog.conf)"
echo "Using intial configuration from $from"
if [[ $from = "maindir" ]]; then
  if [[ -n $conf ]]; then
    cp_from_main_dir "$conf"
  else
    echo "Option cg.inverse.$sim_prog.conf was empty, so I assume $sim_prog needs no conf or you have added it to cg.inverse.filelist."
  fi
elif [[ $from = "laststep" ]]; then
  [[ -n $conf ]] && die "${0##*/}: for initial_configuration '$from' option cg.inverse.$sim_prog.conf is needed!"
  confout="$(csg_get_property --allow-empty cg.inverse.$sim_prog.conf_out)"
  [[ -n $confout ]] && die "${0##*/}: for initial_configuration '$from' option cg.inverse.$sim_prog.confout is needed!"
  #avoid overwriting $confout
  cp_from_last_step --rename "${confout}" "${conf}"
else
  die "${0##*/}: initial_configuration '$from' not implemented"
fi

#convert potential in format for sim_prog
for_all "non-bonded bonded" do_external convert_potential ${sim_prog} '$(csg_get_interaction_property name).pot.cur' '$'"(csg_get_interaction_property inverse.$sim_prog.table)"
