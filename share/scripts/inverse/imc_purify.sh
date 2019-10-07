#!/bin/bash
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
This scripts cleans up the dpot tables for each interaction when using IMC

Usage: ${0##*/}
EOF
   exit 0
fi

name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)
step=$(csg_get_interaction_property step)
kBT=$(csg_get_property cg.inverse.kBT)
bondtype="$(csg_get_interaction_property bondtype)"

echo "purifying dpot for $name"

comment="$(get_table_comment)"
critical csg_resample --in ${name}.dpot.imc --out ${name}.dpot.impure --grid ${min}:${step}:${max} --comment "$comment"
step_nr=$(get_current_step_nr)
scheme=( $(csg_get_interaction_property inverse.do_potential) )
scheme_nr=$(( ( $step_nr - 1 ) % ${#scheme[@]} ))

if [ "${scheme[$scheme_nr]}" = 1 ]; then
  echo "Update potential ${name} : yes"
  do_external table linearop --withflag o ${name}.dpot.impure ${name}.dpot.impure 0 0
  do_external table linearop --withflag i ${name}.dpot.impure ${name}.dpot.impure $kBT 0

  do_external potential shift --type "${bondtype}" ${name}.dpot.impure ${name}.dpot.new
else
  echo "Update potential ${name} : no"
  do_external table linearop ${name}.dpot.impure ${name}.dpot.new 0 0
fi

