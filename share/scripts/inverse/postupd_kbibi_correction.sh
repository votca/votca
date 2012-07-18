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
This script implemtents the post update routine for
the various Kirkwood-Buff corrections

Usage: ${0##*/}
EOF
   exit 0
fi

step_nr=$(get_current_step_nr)
name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)
step=$(csg_get_interaction_property step)

[[ $(csg_get_interaction_property bondtype) = "non-bonded" ]] || die "${0##*/}: kbibi correction only makes sense for non-bonded interactions!"

# always calculate the kbint as there could be cross interaction changes
# needs current rdf and target rdf
if [[ ! -f ${name}.dist.new ]]; then
  do_external rdf $(csg_get_property cg.inverse.program)
fi
if [[ ! -f ${name}.dist.tgt ]]; then
  do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt"
fi
do_external calc kbint ${name}.dist.tgt ${name}.kbint.tgt
if [[ $(csg_get_interaction_property inverse.post_update_options.kbibi.kbint_with_errors) = "yes" ]]; then
  sim_prog="$(csg_get_property cg.inverse.program)"
  rdf_with_errors = $(csg_get_property cg.inverse.${sim_prog}.rdf.with_errors) #filter me away
  [[ ${rdf_with_errors} != "yes" ]] && die "${0##*/}: kb integrals with errors need cg.inverse.${sim_prog}.rdf.with_errors to be yes"
  for f in ${name}_*.dist.block; do
    [[ -f $f ]] || die "${0##*/}: rdf block (${name}_*.dist.block) files not found"
    do_external calc kbint ${f} ${f%.dist.block}.kbint.block
  done
  do_external table average --output ${name}.kbint.new ${name}_*.kbint.block
else
  do_external calc kbint ${name}.dist.new ${name}.kbint.new
fi

kbibi=( $(csg_get_interaction_property inverse.post_update_options.kbibi.do) )
kbibi_nr=$(( ($step_nr - 1 ) % ${#kbibi[@]} ))
if [[ ${kbibi[$kbibi_nr]} = 1 ]]; then
   echo "Apply kbibi correction for interaction ${name}"
   kbibi_type="$(csg_get_interaction_property inverse.post_update_options.kbibi.type)"
   tmpfile=$(critical mktemp ${name}.kbibi.XXX)
   if [[ ${kbibi_type} = "ramp" ]]; then
     do_external kbibi ramp_correction "${name}.kbint.tgt" "${name}.kbint.new" "${tmpfile}"
   elif [[ ${kbibi_type} = "integral" ]]; then
     do_external table integrate --sphere --from left "${name}.dist.tgt" "${name}.dist.tgt.int"
     do_external table integrate --sphere --from left "${name}.dist.new" "${name}.dist.new.int"
     do_external update ibi_pot "${name}.dist.tgt.int" "${name}.dist.new.int" "${name}.pot.cur" "${tmpfile}"
   else
     die "${0##*/}: kbibi type $kbibi_type not implemented yet"
   fi
   comment="$(get_table_comment ${tmpfile})"
   tmpfile2=$(critical mktemp ${name}.kbibi.resample.XXX)
   critical csg_resample --in "${tmpfile}" --out "${tmpfile2}" --grid $min:$step:$max --comment "$comment"
   do_external table add "$1" "${tmpfile2}" "$2"
else
   echo "NO kbibi correction for interaction ${name}"
   do_external postupd dummy "$1" "$2"
fi
