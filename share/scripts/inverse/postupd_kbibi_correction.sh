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
This script implemtents the function update for a single pair
for the Inverse Boltzmann Method

Usage: ${0##*/}
EOF
   exit 0
fi

step_nr=$(get_current_step_nr)
name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)
step=$(csg_get_interaction_property step)

ramp=( $(csg_get_interaction_property inverse.post_update_options.kbibi.do) )
ramp_nr=$(( ($step_nr - 1 ) % ${#ramp[@]} ))
if [ "${ramp[$ramp_nr]}" = 1 ]; then
   echo "Apply kbibi correction for interaction ${name}"
   # needs current rdf and target rdf
   if [[ ! -f ${name}.dist.new ]]; then
     do_external rdf $(csg_get_property cg.inverse.program)
   fi
   if [[ ! -f ${name}.dist.tgt ]]; then
     do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt"
   fi
   tmpfile=$(critical mktemp ${name}.kbibi.XXX)
   do_external kbibi correction ${name}.dist.tgt ${name}.dist.new ${tmpfile}
   tmpfile2=$(critical mktemp ${name}.kbibi.resample.XXX)
   comment="$(get_table_comment ${name}.pressure_correction)"
   critical csg_resample --in ${tmpfile} --out ${tmpfile2} --grid $min:$step:$max --comment "$comment"
   do_external table add "$1" ${tmpfile2} "$2"
else
   echo "NO kbibi correction for interaction ${name}"
   do_external postupd dummy "$1" "$2"
fi
