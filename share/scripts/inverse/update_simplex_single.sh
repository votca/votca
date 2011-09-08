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
- calculates the new property
- compares it to the target property and calculates the target function accordingly

Usage: ${0##*/}
EOF
   exit 0
fi

name=$(csg_get_interaction_property name)
targets=( $(csg_get_interaction_property inverse.simplex.targets) )
weights=( $(csg_get_interaction_property inverse.simplex.weights 1) )
[[ ${#targets[@]} -eq ${#weights[@]} ]] || die "${0##*/}: Number of targets (${#targets[@]}) differ from number of weights (${#weights[@]})"
sim_prog="$(csg_get_property cg.inverse.program)"

sum=0
for ((i=0;i<${#targets[@]};i++)); do
  do_external simplex_target "${targets[$i]}"
  out="${name}.${targets[$i]}.conv"
  [[ -f "${out}" ]] || die "${0##*/}: Could not find '${out}'"
  val="$(<$out)"
  is_num "$val" || die "${0##*/}: Content of $out was not a number"
  x=$(csg_calc ${val} "*" "${weights[$i]}")
  sum=$(csg_calc "$sum" + "$x")
done
echo "$sum" > "${name}.conv"
