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

This script generates a single potential (.pot.new) out a parameter value string (1st argument)

Usage: ${0##*/} parametervalues
EOF
  exit 0
fi

[[ -z $1 ]] && die "${0##*/}: missing argument"
[[ $(echo "$@" | sed -n '$=') -ne 1 ]] && die "${0##*/}: arguments are longer than one line"
n_p="$(echo "$@" | critical awk -F '@' '{print NF}')"

names=( $(csg_get_property cg.non-bonded.name) )
[[ $(( $n_p - 1 )) -ne ${#names[@]} ]] && die "${0##*/}: length of parameter string ($n_p) does not match number of interactions (${#names[@]})"

name="$(csg_get_interaction_property name)"
for ((i=0;i<${#names[@]};i++)); do
  [[ $name = ${names[$i]} ]] && break
done
[[ $name = ${names[$i]} ]] || die "${0##*/}: Could not find interaction $nane in list of all interactions ${names[@]}"

values=( $(echo "$@" | critical awk -v x=$(($i+1)) -F '@' '{print $x}') )
parameters=( $(csg_get_interaction_property inverse.simplex.parameters) )
[[ ${#values[@]} -ne ${#parameters[@]} ]] && die "${0##*/}: Number of values(${#values[@]}) mismatch number of parameters given in xml file (${#parameters[@]})"

para=()
for ((i=0;i<${#values[@]};i++)); do
  para[${#para[@]}]="--var"
  para[${#para[@]}]="${parameters[$i]}=${values[$i]}"
done

fct=$(csg_get_interaction_property inverse.simplex.function)
header=$(csg_get_interaction_property --allow-empty inverse.simplex.functionfile)
min="$(csg_get_interaction_property min)"
step="$(csg_get_interaction_property step)"
max="$(csg_get_interaction_property max)"

if [[ -z $header ]]; then
  do_external table functional --output "${name}.pot.new" "${para[@]}" --fct "${fct}" --grid "${min}:${step}:${max}"
else
  do_external table functional --output "${name}.pot.new" "${para[@]}" --fct "${fct}" --grid "${min}:${step}:${max}" --headerfile "${header}"
fi



