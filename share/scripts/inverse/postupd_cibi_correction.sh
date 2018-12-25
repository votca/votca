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
the integral Kirkwood-Buff corrections described in:
T. E. de Oliveira, P. A. Netz, K. Kremer, C. Junghans, and D. Mukherji,
C-IBI: Targeting cumulative coordination within an iterative protocol
to derive coarse-grained models of (multi-component) complex fluids,
J. Chem. Phys. (in press).

Usage: ${0##*/}
EOF
   exit 0
fi

step_nr=$(get_current_step_nr)
name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)
step=$(csg_get_interaction_property step)

[[ $(csg_get_interaction_property bondtype) = "non-bonded" ]] || die "${0##*/}: cibi correction only makes sense for non-bonded interactions!"

# always calculate the kbint as there could be cross interaction changes
# needs current rdf and target rdf
if [[ ! -f ${name}.dist.new ]]; then
  do_external rdf $(csg_get_property cg.inverse.program)
fi
if [[ ! -f ${name}.dist.tgt ]]; then
  do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt"
fi
do_external calc kbint ${name}.dist.tgt ${name}.kbint.tgt
if [[ $(csg_get_interaction_property inverse.post_update_options.cibi.kbint_with_errors) = "yes" ]]; then
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

cibi=( $(csg_get_interaction_property inverse.post_update_options.cibi.do) )
cibi_nr=$(( ($step_nr - 1 ) % ${#cibi[@]} ))
if [[ ${cibi[$cibi_nr]} = 1 ]]; then
   echo "Apply cibi correction for interaction ${name}"
   tmpfile=$(critical mktemp ${name}.cibi.XXX)
   do_external table integrate --sphere --from left "${name}.dist.tgt" "${name}.dist.tgt.int"
   do_external table integrate --sphere --from left "${name}.dist.new" "${name}.dist.new.int"
   kBT="$(csg_get_property cg.inverse.kBT)"
   is_num "${kBT}" || die "${0##*/}: cg.inverse.kBT should be a number, but found '$kBT'"
   do_external update ibi_pot "${name}.dist.tgt.int" "${name}.dist.new.int" "${name}.pot.cur" "${tmpfile}" "${kBT}"
   comment="$(get_table_comment ${tmpfile})"
   tmpfile2=$(critical mktemp ${name}.cibi.resample.XXX)
   critical csg_resample --in "${tmpfile}" --out "${tmpfile2}" --grid $min:$step:$max --comment "$comment"
   do_external table add "$1" "${tmpfile2}" "$2"
else
   echo "No cibi correction for interaction ${name}"
   do_external postupd dummy "$1" "$2"
fi
