#! /bin/bash
#
# Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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
This script does a Gauss-Newton update of non-bonded interactions.

It expects to find jacobian.npz already generated

Usage: ${0##*/}  [--help]
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"

# pressure constraint
pressure_constraint=$(csg_get_property cg.inverse.gauss_newton.pressure_constraint)
if is_num "${pressure_constraint}"; then
  p_file="${name}.pressure"
  do_external pressure "$sim_prog" "$p_file"
  p_now="$(sed -n 's/^Pressure=\(.*\)/\1/p' "$p_file")" || die "${0##*/}: sed of Pressure failed"
  [[ -z $p_now ]] && die "${0##*/}: Could not get pressure from simulation"
  echo "New pressure $p_now, target pressure $pressure_constraint"
  pressure_constraint_flag="--pressure-constraint ,$pressure_constraint,$p_now"
fi

# Kirkwood-Buff integral constraint
if [[ $(csg_get_property cg.inverse.gauss_newton.kirkwood_buff_constraint) == 'true' ]]; then
    kirkwood_buff_constraint_flag="--kirkwood-buff-constraint"
fi

# potential energy constraint
potential_energy_constraint=$(csg_get_property cg.inverse.gauss_newton.potential_energy_constraint)
if is_num "${potential_energy_constraint}"; then
  PE_file="${name}.potential_energy"
  do_external potential_energy "$sim_prog" "$PE_file"
  PE_now="$(sed -n 's/^Potential=\(.*\)/\1/p' "$p_file")" || die "${0##*/}: sed of Potential failed"
  [[ -z $PE_now ]] && die "${0##*/}: Could not get potential energy from simulation"
  echo "New potential energy $PE_now, target potential energy $potential_energy_constraint"
  pressure_constraint_flag="--potential-energy-constraint ,$potential_energy_constraint,$PE_now"
fi

# Gauss-Newton residual weighting
residual_weighting_flag="--residual-weighting '$(csg_get_property cg.inverse.gauss_newton.residual_weighting)'"

# topology for molecular conections and volume
topol=$(csg_get_property cg.inverse.topol_xml)
[[ -f $topol ]] || die "${0##*/}: topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist"

# volume
volume=$(critical csg_dump --top "$topol" | grep 'Volume' | awk '{print $2}')
([[ -n "$volume" ]] && is_num "$volume") || die "could not determine the volume from file ${topol}"

# verbose
verbose=$(csg_get_property cg.inverse.verbose)
step_nr=$(get_current_step_nr)
[[ "${verbose}" == 'true' ]] && verbose_flag="--verbose"
[[ "${verbose}" == 'step0+1' ]] && [[ $step_nr == '0' || $step_nr == '1' ]] && verbose_flag="--verbose"

# which interactions to update and target
export step_nr=$(get_current_step_nr)  # needs to be exported to work in for_all
do_potential_list="$(for_all "non-bonded" 'scheme=( $(csg_get_interaction_property inverse.do_potential) ); echo -n $(csg_get_interaction_property name),${scheme[$(( ($step_nr - 1 ) % ${#scheme[@]} ))]}:')"
tgt_dist_list="$(for_all "non-bonded" 'tgt_dist=( $(csg_get_interaction_property inverse.is_target_distribution) ); echo -n $(csg_get_interaction_property name),${tgt_dist}:')"
upd_pots_flag="--upd-pots $do_potential_list"
tgt_dists_flag="--tgt-dists $tgt_dist_list"

# weather to make dU = 0 one point before the cut off
[[ "$(csg_get_property cg.inverse.gauss_newton.flatten_at_cut_off)" == 'true' ]] && flatten_at_cut_off_flag="--flatten-at-cut-off"

# Some arguments (cut_off, kBT) will be read directly from the settings.xml. They do not have a default in csg_defaults.xml.
# Others (closure, ...) could also be read from the settings file, but this bash script handles the defaults.
do_external update gauss_newton_py gauss-newton \
  ${verbose_flag-} \
  --jacobian "jacobian.npz" \
  --volume "$volume" \
  --topol "$topol" \
  --options "$CSGXMLFILE" \
  --g-tgt-ext "dist.tgt" \
  --g-cur-ext "dist.new" \
  ${pressure_constraint_flag-} \
  ${kirkwood_buff_constraint_flag-} \
  ${residual_weighting_flag-} \
  ${upd_pots_flag-} \
  ${tgt_dists_flag-} \
  ${flatten_at_cut_off_flag-} \
  --out "dpot.pure_gn"

# overwrite with zeros if do_potential=0
do_potential_zero_overwrite() {
  step_nr=$(get_current_step_nr)
  scheme=( $(csg_get_interaction_property inverse.do_potential) )
  scheme_nr=$(( ( $step_nr - 1 ) % ${#scheme[@]} ))
  name=$(csg_get_interaction_property name)
  if [[ ${scheme[$scheme_nr]} == 0 ]]; then
    echo "Update potential ${name} : no"
    min=$(csg_get_interaction_property min)
    max=$(csg_get_interaction_property max)
    step=$(csg_get_interaction_property step)
    critical rm "${name}.dpot.new"
    do_external table dummy "${min}:${step}:${max}" "${name}.dpot.new"
  fi
}
export -f do_potential_zero_overwrite
for_all "non-bonded" do_potential_zero_overwrite
