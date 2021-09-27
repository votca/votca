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
This script generates the initial state file and puts all in-file together

Usage: ${0##*/} outputfile
EOF
   exit 0
fi

[[ -z $1 ]] && die "${0##*/}: Missing argument"

names="$(csg_get_interaction_property --all name)"
parameters="$(csg_get_interaction_property --all inverse.optimizer.parameters)"
what=$(has_duplicate "${parameters[@]}") && die "${0##*/}: the parameter $what appears twice"
otype="$(csg_get_property cg.inverse.optimizer.type)"

liste=()
for name in $names; do
  input="$name.${otype}.in"
  [[ -f $input ]] || die "${0##*/}: Could not find $input"
  [[ -n $(sed -n '/^[#@]/p' "$input") ]] && die "${0##*/}: $input has commentlines"
  liste[${#liste[@]}]="$input"
done

get_table_comment | sed 's/^/#/' > "$1"
echo "#Interactions: ${names}" >> "$1"
echo "#Method: ${otype}" >> "$1"
echo "#State = Initialization" >> $1
echo "#Format $parameters conv flag" >> "$1"
#added conv=0 and flag=pending to all lines
critical paste "${liste[@]}" | critical sed -e 's/$/ 0 pending/' >> "$1"

#cma converts initial values to first population
if [[ $otype = cma ]]; then
  tmpfile="$(critical mktemp "$1.XXX")"
  eps="$(csg_get_property cg.inverse.optimizer.cma.eps)"
  critical mv "$1" "$tmpfile"
  do_external cma precede_state --eps "$eps" "$tmpfile" "$1"
fi

