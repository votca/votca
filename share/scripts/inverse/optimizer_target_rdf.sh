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
do_external rdf ${sim_prog}
name="$(csg_get_interaction_property name)"
[[ -f ${name}.dist.new ]] || die "${0##*/}: Could not calculate ${name}.dist.new"
target="$(csg_get_interaction_property inverse.optimizer.rdf.target)"
do_external resample target ${target} ${name}.dist.tgt
weight=$(csg_get_interaction_property --allow-empty inverse.optimizer.rdf.weight)
weightfile=$(csg_get_interaction_property --allow-empty inverse.optimizer.rdf.weightfile)
if [[ -n $weight ]]; then
  min=$(csg_get_interaction_property min)
  step=$(csg_get_interaction_property step)
  max=$(csg_get_interaction_property max)
  if [[ -z ${weightfile} ]]; then
    do_external table functional --fct "$weight" --grid "${min}:${step}:${max}" "${name}.rdf.weight"
  else
    [[ -f ${weightfile} ]] || die "${0##*/}: Could not find function weightfile '$weightfile'. did you forget to add it to cg.inverse.filelist?"
    do_external table functional --fct "$weight" --grid "${min}:${step}:${max}" --headerfile "${weightfile}" "${name}.rdf.weight"
  fi
  do_external table combine --op x "${name}.dist.tgt" "${name}.rdf.weight" "${name}.dist.tgt.weighted"
  do_external table combine --op x "${name}.dist.new" "${name}.rdf.weight" "${name}.dist.new.weighted"
  do_external table combine --sum --op d "${name}.dist.tgt.weighted" "${name}.dist.new.weighted" > "${name}.rdf.conv"
else
  do_external table combine --sum --op d "${name}.dist.tgt" "${name}.dist.new" > "${name}.rdf.conv"
fi

