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
This script does a multistate Gauss-Newton update of non-bonded interactions.

It expects to find jacobian.npz per state already generated

Usage: ${0##*/}  [--help]
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
state_names="$(csg_get_property cg.inverse.multistate.state_names)"
state_names_arr=( $(csg_get_property cg.inverse.multistate.state_names) )
state_kBTs_arr=( $(csg_get_property cg.inverse.multistate.state_kBTs) )

# verbose
verbose=$(csg_get_property cg.inverse.verbose)
step_nr=$(get_current_step_nr)
[[ "${verbose}" == 'true' ]] && verbose_flag="--verbose"
[[ "${verbose}" == 'step0+1' ]] && [[ $step_nr == '0' || $step_nr == '1' ]] && verbose_flag="--verbose"

# no constraints
pressure_constraint="$(csg_get_property cg.inverse.gauss_newton.pressure_constraint)"
kbi_constraint="$(csg_get_property cg.inverse.gauss_newton.kirkwood_buff_constraint)"
potential_energy_constraint=$(csg_get_property cg.inverse.gauss_newton.potential_energy_constraint)
if is_num "${pressure_constraint}"; then
  die "No pressure constraint implemented for multistate"
elif [[ $kbi_constraint == 'true' ]]; then
  die "No KBI constraint implemented for multistate"
elif is_num "${potential_energy_constraint}"; then
  die "No PE constraint implemented for multistate"
fi

# Gauss-Newton residual weighting
residual_weighting_flag="--residual-weighting '$(csg_get_property cg.inverse.gauss_newton.residual_weighting)'"

# topology, volume, and jacobian
jacobian_flag="--jacobian"
for state in $state_names; do
  topol_state="${state}/$(csg_get_property cg.inverse.topol_xml)"
  [[ -f $topol_state ]] || die "${0##*/}: topol file '$topol_state' not found, possibly you have to add it to cg.inverse.filelist"
  volume_state=$(critical csg_dump --top "$topol_state" | grep 'Volume' | awk '{print $2}')
  ([[ -n "$volume_state" ]] && is_num "$volume_state") || die "could not determine the volume from file ${topol_state}"
  # append
  topol="$topol $topol_state"
  volume="$volume $volume_state"
  jacobian_flag="$jacobian_flag ${state}/jacobian.npz"
done

# which interactions to update and target
export step_nr=$(get_current_step_nr)  # needs to be exported to work in for_all
do_potential_list="$(for_all "non-bonded" 'scheme=( $(csg_get_interaction_property inverse.do_potential) ); echo -n $(csg_get_interaction_property name),${scheme[$(( ($step_nr - 1 ) % ${#scheme[@]} ))]}:')"
tgt_dist_list="$(for_all "non-bonded" 'tgt_dist=( $(csg_get_interaction_property inverse.is_target_distribution) ); echo -n $(csg_get_interaction_property name),${tgt_dist}:')"
upd_pots_flag="--upd-pots $do_potential_list"
tgt_dists_flag="--tgt-dists $tgt_dist_list"

# weather to make dU = 0 one point before the cut off
[[ "$(csg_get_property cg.inverse.gauss_newton.flatten_at_cut_off)" == 'true' ]] && flatten_at_cut_off_flag="--flatten-at-cut-off"

# Some arguments (cut_off, kBT, state_kBTs, ..) will be read directly from the settings.xml. They do not have a default in csg_defaults.xml.
# Others (upd_pots, tgt_dists, ...) could also be read from the settings file, but this bash script handles the defaults.
do_external update gauss_newton_py gauss-newton \
  ${verbose_flag-} \
  --multistate \
  ${jacobian_flag-} \
  --volume $volume \
  --topol $topol \
  --options "$CSGXMLFILE" \
  --g-tgt-ext "dist.tgt" \
  --g-cur-ext "dist.new" \
  ${residual_weighting_flag-} \
  ${tgt_dcdh_flag-} \
  ${upd_pots_flag-} \
  ${tgt_dists_flag-} \
  ${flatten_at_cut_off_flag-} \
  --out "dpot.new"

# resample potentials. This is needed because non-bonded.max is usually larger than iie.cut-off and the former should define the table
for_all "non-bonded" 'csg_resample --in $(csg_get_interaction_property name).dpot.new --out $(csg_get_interaction_property name).dpot.new --grid $(csg_get_interaction_property min):$(csg_get_interaction_property step):$(csg_get_interaction_property max) --comment "adapted to grid in update_iie.sh"'
# csg_resample alone will not do the job
for_all "non-bonded" 'do_external table extrapolate --function constant --region right $(csg_get_interaction_property name).dpot.new $(csg_get_interaction_property name).dpot.new'

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
