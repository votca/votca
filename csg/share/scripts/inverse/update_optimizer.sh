#! /usr/bin/env bash
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
- performs optimizer algorithm if no pending parameter sets present
- continues with next parameter set in table if otherwise

Usage: ${0##*/}
EOF
   exit 0
fi

otype="$(csg_get_property cg.inverse.optimizer.type)"

#generate mappings to calculate reference rdf, which depend on all of them
for_all "non-bonded bonded" do_external optimizer state_to_mapping "${otype}.state.cur"

for_all "non-bonded bonded" do_external update optimizer_single

names="$(csg_get_interaction_property --all name)"
conv=0
for name in ${names}; do
  [[ -f ${name}.conv ]] || die "${0##*/}: could not find '${name}.conv'"
  x=$(<${name}.conv)
  is_num "$x" || die "${0##*/}: content of '${name}.conv' was not a number"
  conv=$(csg_calc "$conv" + "$x")
done

[[ -f ${otype}.state.cur ]] || die "${0##*/}: Could not find ${otype}.state.cur"
active="$(critical sed -n '/0 active$/{=;q}' "${otype}.state.cur")"
[[ -z $active ]] && die "${0##*/}: not could find an active simulation in ${otype}.state.cur"
is_int "$active" || die "${0##*/}: Strange - $active should be a number"

#check if there are still pending simulations
pending="$(critical sed -n '/pending/p' "${otype}.state.cur")"
if [[ -z $pending ]]; then
  #simplex needs to know which one was the last try guess and it need $conv, cma does not care
  critical sed "${active}s/0 active$/ $conv try/" "${otype}.state.cur" > "${otype}.state.try"
  do_external ${otype} precede_state "${otype}.state.try" "${otype}.state.done"
  state="$(critical sed -n 's/#State = \(.*\)$/\1/p' ${otype}.state.done)"
  [[ ${otype} = simplex ]] && msg "Simplex state changed to $state"
else
  pending="$(echo "$pending" | critical sed -n '$=')"
  msg "There are still $pending simulations to be performed before the next ${otype} state change"
  critical sed "${active}s/0 active$/ $conv complete/" "${otype}.state.cur" > "${otype}.state.done"
  if [[ ${otype} = cma ]]; then
    critical cp "${otype}.internal_state.cur" "${otype}.internal_state.new"
  fi
fi
do_external optimizer state_to_potentials "${otype}.state.done" "${otype}.state.new"
