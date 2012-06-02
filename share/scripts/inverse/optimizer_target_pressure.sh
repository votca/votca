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

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
Calculates the difference current and target pressure

Usage: ${0##*/}
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
name="$(csg_get_interaction_property name)"

p_file="${name}.pressure"
do_external pressure "$sim_prog" "$p_file" 
p_now="$(sed -n 's/^Pressure=\(.*\)/\1/p' "$p_file")" || die "${0##*/}: sed of Pressure failed"
[[ -z $p_now ]] && die "${0##*/}: Could not get pressure from simulation"
echo "New pressure $p_now"

if [[ $p_now = nan ]]; then
  p_undef=$(csg_get_interaction_property --allow-empty cg.inverse.optimizer.pressure.undef)
  [[ -z $p_undef ]] && die "${0##*/}: Pressure was nan, check your simulation (this usually means system has blow up -> use pre simulation)"
  p_now="${p_undef}"
fi

p_target="$(csg_get_interaction_property inverse.p_target)"

critical awk -v x="$p_now" -v y="$p_target" 'BEGIN{print sqrt((x-y)^2)}' > "${name}.pressure.conv"

