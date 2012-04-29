#! /bin/bash
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

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script implements the initialization for every step in a generic way

Usage: ${0##*/}
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
otype="$(csg_get_property cg.inverse.optimizer.type)"

#get new pot from last step and make it current potential
for_all "non-bonded bonded" 'cp_from_last_step --rename $(csg_get_interaction_property name).pot.new $(csg_get_interaction_property name).pot.cur'

cp_from_last_step --rename "${otype}.state.new" "${otype}.state.cur"
if [[ $otype = cma ]]; then
  cp_from_last_step --rename "${otype}.internal_state.new" "${otype}.internal_state.cur"
fi

#initialize sim_prog
do_external initstep_generic $sim_prog
