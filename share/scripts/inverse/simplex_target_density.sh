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
Calculated the difference between rdf

Usage: ${0##*/}
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
name="$(csg_get_interaction_property name)"
mol="$(csg_get_interaction_property inverse.simplex.density.molname "*")"
axis="$(csg_get_interaction_property inverse.simplex.density.axis "x")"
step="$(csg_get_interaction_property inverse.simplex.density.step)"
opts=( "--molname" "$mol" "--axis" "$axis" "--step" "$step" )
do_external density ${sim_prog} "${name}.density.new" "${opts[@]}"

[[ -f ${name}.density.new ]] || die "${0##*/}: Could not calculate ${name}.density.new"
target="$(csg_get_interaction_property inverse.simplex.density.target)"
main_dir=$(get_main_dir)
min="$(csg_get_interaction_property inverse.simplex.density.min)"
max="$(csg_get_interaction_property inverse.simplex.density.max)"
step="$(csg_get_interaction_property inverse.simplex.density.step)"

[[ -f ${main_dir}/$target ]] || die "${0##*/}: Not find density target '$target' in maindir"
critical csg_resample --in "${main_dir}/${target}" --out "${name}.density.tgt" --grid "$min:$step:$max"
t1=$(critical mktemp "${name}.density.new.cut.XXXX")
critical csg_resample --in ${name}.density.new --out "$t1" --grid "$min:$step:$max"
do_external table combine --sum --no-flags --op d "${t1}" "${name}.density.tgt" > "${name}.density.conv"

