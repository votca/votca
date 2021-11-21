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

This script
- reads sinple interaction optimizer infile
- checks if the number of values are enough

Usage: ${0##*/} N

where N is the total number of parameters
EOF
  exit 0
fi

[[ -z $1 ]] && die "${0##*/}: missing argument"
is_int "$1" || die "${0##*/}: 1st argument should be a int"
total_parameter=$(( $1 + 1 ))
parameters=( $(csg_get_interaction_property inverse.optimizer.parameters) )
name=$(csg_get_interaction_property name)
otype="$(csg_get_property cg.inverse.optimizer.type)"

[[ ${#parameters[@]} -gt 0 ]] || die "${0##*/}: property inverse.optimizer.parameters for interaction $name should contain more than 1 parameter"

main_dir=$(get_main_dir)
input="$name.${otype}.in"

if [[ -f $main_dir/$input ]]; then
  msg "Using ${otype} parameter input ${input} for interaction ${name}"
  cp_from_main_dir $input
  critical sed -i '/^[#@]/d' "$input"
  n=( $(critical awk '{print NF}' "$input") )
  if [[ ${otype} = simplex ]]; then
    [[ ${#n[@]} -eq $total_parameter ]] || die "${0##*/}: Simplex needs exactly $total_parameter parameter sets (lines) in $input, but got ${#n[@]} lines"
  elif [[ ${otype} = cma ]]; then
    [[ ${#n[@]} -eq 1 ]] || die "${0##*/}: cma needs exactly 1 parameter sets (lines) in $input"
  else
    die "${0##*/} Optimizer method '${otype}' not implemented, yet!"
  fi
  for ((i=0;i<${#n[@]};i++)); do
    [[ ${n[$i]} -eq ${#parameters[@]} ]] || die "${0##*/}: Number of paramenter in line $(($i+1)) of $input is ${n[$i]} and so unequal to the number given in xml file (${#parameters[@]})"
  done
else
  die "No input file $input found"
fi
