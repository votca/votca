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

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script implemtents the function update for a single pair
for the Inverse Boltzmann Method

Usage: ${0##*/}
EOF
   exit 0
fi

step_nr=$(get_current_step_nr)
sim_prog="$(csg_get_property cg.inverse.program)"
scheme=( $(csg_get_interaction_property inverse.do_potential) )
scheme_nr=$(( ($step_nr - 1 ) % ${#scheme[@]} ))
name=$(csg_get_interaction_property name)
bondtype="$(csg_get_interaction_property bondtype)"
extrap_near_core=$(csg_get_property cg.inverse.hncgn.extrap_near_core)
fix_near_cut_off=$(csg_get_property cg.inverse.hncgn.fix_near_cut_off)

pressure_constraint=$(csg_get_property cg.inverse.hncgn.pressure_constraint)
if is_num ${pressure_constraint}; then
    p_file="${name}.pressure"
    do_external pressure "$sim_prog" "$p_file"
    p_now="$(sed -n 's/^Pressure=\(.*\)/\1/p' "$p_file")" || die "${0##*/}: sed of Pressure failed"
    [[ -z $p_now ]] && die "${0##*/}: Could not get pressure from simulation"
    echo "New pressure $p_now, target pressure $p_target"
fi

if [ "${scheme[$scheme_nr]}" = 1 ]; then
    echo "Update potential ${name} : yes"
    #update hncgn
    do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt"
    kBT="$(csg_get_property cg.inverse.kBT)"
    density="$(csg_get_property cg.inverse.hncgn.density)"
    cut_off="$(csg_get_property cg.inverse.hncgn.cut_off)"
    is_num "${kBT}" || die "${0##*/}: cg.inverse.kBT should be a number, but found '$kBT'"
    if is_num ${pressure_constraint}; then
        do_external update hncgn_pot ${name}.dist.tgt ${name}.dist.new ${name}.pot.cur ${name}.dpot.pure_hncgn "${kBT}" "${density}" "${cut_off}" \
        "--pressure-constraint=${pressure_constraint},${p_now}" "--extrap-near-core=${extrap_near_core}" "--fix-near-cut-off=${fix_near_cut_off}"
    else
        do_external update hncgn_pot ${name}.dist.tgt ${name}.dist.new ${name}.pot.cur ${name}.dpot.pure_hncgn "${kBT}" "${density}" "${cut_off}" \
        "--extrap-near-core=${extrap_near_core}" "--fix-near-cut-off=${fix_near_cut_off}"
    fi
   do_external potential shift --type "${bondtype}" ${name}.dpot.pure_hncgn ${name}.dpot.new
else
   echo "Update potential ${name} : no"
   min=$(csg_get_interaction_property min)
   step=$(csg_get_interaction_property step)
   max=$(csg_get_interaction_property max)
   do_external table dummy ${min}:${step}:${max} ${name}.dpot.new
fi
