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
This script:
- implements the update function for each non-bonded interaction
- performs downhill simplex algorithm if no pending parameter sets present
- continues with next parameter set in table if otherwise

Usage: ${0##*/}
EOF
   exit 0
fi

for_all non-bonded do_external update simplex_single

names="$(csg_get_property cg.non-bonded.name)"
conv=0
for name in ${names}; do
  [[ -f ${name}.conv ]] || die "${0##*/}: could not find '${name}.conv'"
  x=$(<${name}.conv)
  is_num "$x" || die "${0##*/}: content of '${name}.conv' was not a number"
  conv=$(csg_calc "$conv" + "$x")
done

[[ -f simplex.table.cur ]] || die "${0##*/}: Could not find simplex.table.cur"
active="$(critical sed -n '/active$/{=;q}' "simplex.table.cur")"
[[ -z $active ]] && die "${0##*/}: not could find an active simulation in simplex.table.cur"
is_int "$active" || die "${0##*/}: Strange - $active should be a number"

#check if there are still pending simulations
pending="$(critical sed -n '/pending/p' "simplex.table.done")"
if [[ -z $pending ]]; then
  #simplex needs to know which one was the try guess
  critical sed "${active}s/[^#]*$/$conv try/" "simplex.table.cur" > "simplex.table.try"
  do_external simplex precede_state "simplex.state.cur" "simplex.table.try" "simplex.state.new"
  do_external simplex update_table "simplex.state.new" "simplex.table.try" "simplex.table.done"
else
  critical sed "${active}s/[^#]*$/$conv complete/" "simplex.table.cur" > "simplex.table.done"
  pending="$(echo "$pending" | critical sed -n '$=')"
  msg "There are still $pending simulations to be performed before the next simplex state change"
  critical cp "simplex.state.cur" "simplex.state.new"
fi
do_external simplex table_to_potentials "simplex.table.done" "simplex.table.new"
