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

This script generates a mapping for the reference mapping from the parameters of the active in input state using the mapping template

Usage: ${0##*/} input
EOF
  exit 0
fi

[[ -z $1 ]] && die "${0##*/}: missing argument"
input="$1"
[[ -f $input ]] || die "${0##*/}: Could not read $input"

name="$(csg_get_interaction_property name)"
if [[ $(csg_get_interaction_property inverse.optimizer.mapping.change) = no ]]; then
  echo "Interaction $name does not change the reference mapping"
  exit 0
fi

line="$(critical sed -n '/active$/{=;q}' "$input")"
[[ -z $line ]] && die "${0##*/}: not could find a active line in $input"
is_int "$line" || die "${0##*/}: Strange - $line should be a number"
values=( $(sed -n "${line}p" "$input") )

parameters=( $(csg_get_interaction_property --all inverse.optimizer.parameters) )
[[ $(( ${#values[@]} - 2 )) -ne ${#parameters[@]} ]] && die "${0##*/}: length of parameter string (${#values[@]}) does not match number of interactions (${#parameters[@]})"
what=$(has_duplicate "${parameters[@]}") && die "${0##*/}: the parameter $what appears twice"

template="$(csg_get_interaction_property inverse.optimizer.mapping.template)"
output="$(csg_get_interaction_property inverse.optimizer.mapping.output)"
cp_from_main_dir --rename "${template}" "${output}"

for ((i=0;i<${#parameters[@]};i++)); do
  echo "Replacing @${parameters[$i]}@ with ${values[$i]} in $output"
  critical sed -i -e "s/@${parameters[$i]}@/${values[$i]}/g" "${output}"
done
