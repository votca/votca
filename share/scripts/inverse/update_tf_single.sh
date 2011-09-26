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
if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script implemtents the function update of a single interaction 
for the thermodynamics force iteration

Usage: ${0##*/}
EOF
   exit 0
fi

step_nr=$(get_current_step_nr)
scheme=( $(csg_get_interaction_property inverse.do_potential) )
scheme_nr=$(( ($step_nr - 1 ) % ${#scheme[@]} ))
name=$(csg_get_interaction_property name)

sim_prog="$(csg_get_property cg.inverse.program)"

mol="$(csg_get_interaction_property tf.molname)"
adress_type=$(get_simulation_setting adress_type)
step=$(csg_get_interaction_property step)
opts=()
if [[ $adress_type = "sphere" ]]; then
  echo "Adress type: $adress_type"
  max=$(csg_get_interaction_property tf.spline_end)
  adressc="$(get_simulation_setting adress_reference_coords "0 0 0")"
  ref="$(echo "$adressc" | awk '{if (NF<3) exit 1; printf "[%s,%s,%s]",$1,$2,$3;}')" || die "${0##*/}: we need three numbers in adress_reference_coords, but got '$adressc'"
  axis="r"
  opts=( "--rmax" "$max" "--ref" "$ref" )
else
  echo "Adress type: $adress_type"
  axis="x"
fi
opts=( "${opts[@]}" "--molname" "$mol" "--axis" "$axis" "--step" "$step" )
do_external density $sim_prog "$name.dist.new" "${opts[@]}"
if [ "${scheme[$scheme_nr]}" = 1 ]; then
   echo "Update tf ${name} : yes"
    do_external calc thermforce ${name}.dist.new ${name}.dpot.new
else
   echo "Update tf ${name} : no"
   min=$(csg_get_interaction_property min)
   step=$(csg_get_interaction_property step)
   max=$(csg_get_interaction_property max)
   do_external table dummy ${min}:${step}:${max} ${name}.dpot.new
fi
