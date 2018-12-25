#! /bin/bash
#
# Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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
the ramp Kirkwood-Buff corrections as described in:
P. Ganguly, D. Mukherji, C. Junghans, N. F. A. van der Vegt,
Kirkwood-Buff coarse-grained force fields for aqueous solutions,
J. Chem. Theo. Comp., 8, 1802 (2012), doi:10.1021/ct3000958

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
  rdf_with_errors=$(csg_get_property cg.inverse.$sim_prog.rdf.with_errors)
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
   tmpfile=$(critical mktemp ${name}.kbibi.XXX)
   kBT="$(csg_get_property cg.inverse.kBT)"
   is_num "${kBT}" || die "${0##*/}: cg.inverse.kBT should be a number, but found '$kBT'"
   int_start=$(csg_get_interaction_property inverse.post_update_options.kbibi.start);
   is_num "${int_start}" || die "${0##*/}: interaction property 'inverse.post_update_options.kbibi.start', should be a number, but found '${int_start}'"
   csg_calc "${int_start}" "<" "${min}" && die "${0##*/}: 'inverse.post_update_options.kbibi.start'(${int_start}) is smaller than min (${min}) for interaction '$name'"
   int_stop=$(csg_get_interaction_property inverse.post_update_options.kbibi.stop);
   is_num "${int_stop}" || die "${0##*/}: interaction property 'inverse.post_update_options.kbibi.stop', should be a number, but found '${int_stop}'"
   csg_calc "${int_stop}" ">" "${max}" && die "${0##*/}: 'inverse.post_update_options.kbibi.stop'(${int_stop}) is bigger than max (${max}) for interaction '$name'"
   ramp_factor=$(csg_get_interaction_property inverse.post_update_options.kbibi.factor);
   is_num "${ramp_factor}" || die "${0##*/}: interaction property 'inverse.post_update_options.kbibi.factor', should be a number, but found '${ramp_factor}'"
   r_ramp=$(csg_get_interaction_property --allow-empty inverse.post_update_options.kbibi.r_ramp)
   do_external kbibi ramp_correction "${name}.kbint.tgt" "${name}.kbint.new" "${tmpfile}" "${kBT}" "$min:$step:${r_ramp:-$max}" "${int_start}:${int_stop}" "${ramp_factor}"
   comment="$(get_table_comment ${tmpfile})"
   tmpfile2=$(critical mktemp ${name}.kbibi.resample.XXX)
   critical csg_resample --in "${tmpfile}" --out "${tmpfile2}" --grid $min:$step:$max --comment "$comment"
   do_external table add "$1" "${tmpfile2}" "$2"
else
   echo "No kbibi correction for interaction ${name}"
   do_external postupd dummy "$1" "$2"
fi
