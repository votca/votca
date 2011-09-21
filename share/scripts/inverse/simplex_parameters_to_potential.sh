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

name=$(csg_get_interaction_property name)
parameters=( $(csg_get_property cg.non-bonded.inverse.simplex.parameters) )
[[ $(( $# - 2 )) -ne ${#parameters[@]} ]] && die "${0##*/}: length of parameter string ($#) does not match number of interactions (${#parameters[@]})"
what=$(has_duplicate "${parameters[@]}") && die "${0##*/}: the parameter $what appears twice"

para=()
for ((i=1;i<=${#parameters[@]};i++)); do
  para[${#para[@]}]="--var"
  para[${#para[@]}]="${parameters[$i-1]}=${!i}"
done

fct=$(csg_get_interaction_property inverse.simplex.function)
header=$(csg_get_interaction_property --allow-empty inverse.simplex.functionfile)
min="$(csg_get_interaction_property min)"
step="$(csg_get_interaction_property step)"
max="$(csg_get_interaction_property max)"

if [[ -z $header ]]; then
  do_external table functional "${para[@]}" --fct "${fct}" --grid "${min}:${step}:${max}" "${name}.pot.new"
else
  do_external table functional "${para[@]}" --fct "${fct}" --grid "${min}:${step}:${max}" --headerfile "${header}" "${name}.pot.new"
fi
