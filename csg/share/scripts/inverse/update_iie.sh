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
This script implements the potential update for the iterative integral equation
methods

Usage: ${0##*/}  [--help]
EOF
   exit 0
fi

iie_method="$(csg_get_property cg.inverse.iie.method)"
sim_prog="$(csg_get_property cg.inverse.program)"

# newton_mod not supported any more
if [[ $iie_method == 'newton-mod' ]]; then
  die "Newton-mod method is not supported any more. Use Newton's method!"
fi

# pressure constraint
if [[ $iie_method == 'gauss-newton' ]]; then
  pressure_constraint=$(csg_get_property cg.inverse.iie.pressure_constraint)
  if is_num "${pressure_constraint}"; then
    p_file="${name}.pressure"
    do_external pressure "$sim_prog" "$p_file"
    p_now="$(sed -n 's/^Pressure=\(.*\)/\1/p' "$p_file")" || die "${0##*/}: sed of Pressure failed"
    [[ -z $p_now ]] && die "${0##*/}: Could not get pressure from simulation"
    echo "New pressure $p_now, target pressure $pressure_constraint"
    pressure_constraint_flag="--pressure-constraint ,$pressure_constraint,$p_now"
  fi
fi

# Kirkwood-Buff integral constraint
if [[ $iie_method == 'gauss-newton' && $(csg_get_property cg.inverse.iie.kirkwood_buff_constraint) == 'true' ]]; then
    kirkwood_buff_constraint_flag="--kirkwood-buff-constraint"
fi

# Gauss-Newton residual weighting
if [[ $iie_method == 'gauss-newton' ]]; then
  residual_weighting_flag="--residual-weighting '$(csg_get_property cg.inverse.iie.residual_weighting)'"
fi

# target dc/dh
if [[ $(csg_get_property cg.inverse.iie.tgt_dcdh) == 'true' ]]; then
  tgt_dcdh_flag="--tgt-dcdh $(get_main_dir)/dcdh.npz"
else
  g_intra_flag="--g-cur-intra-ext dist-intra.new"
fi

# topology for molecular conections and volume
topol=$(csg_get_property --allow-empty cg.inverse.iie.topol)
[[ -z $topol ]] && topol=$(csg_get_property cg.inverse.$sim_prog.topol)
[[ -f $topol ]] || die "${0##*/}: topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist"

# volume
volume=$(critical csg_dump --top "$topol" | grep 'Volume' | awk '{print $2}')
([[ -n "$volume" ]] && is_num "$volume") || die "could not determine the volume from file ${topol}"

# verbose
verbose=$(csg_get_property cg.inverse.iie.verbose)
step_nr=$(get_current_step_nr)
[[ "${verbose}" == 'true' ]] && verbose_flag="--verbose"
[[ "${verbose}" == 'step0+1' ]] && [[ $step_nr == '0' || $step_nr == '1' ]] && verbose_flag="--verbose"

# which interactions to update and target
if [[ $iie_method == 'gauss-newton' ]]; then
  export step_nr=$(get_current_step_nr)  # needs to be exported to work in for_all
  do_potential_list="$(for_all "non-bonded" 'scheme=( $(csg_get_interaction_property inverse.do_potential) ); echo -n $(csg_get_interaction_property name),${scheme[$(( ($step_nr - 1 ) % ${#scheme[@]} ))]}:')"
  tgt_dist_list="$(for_all "non-bonded" 'tgt_dist=( $(csg_get_interaction_property inverse.is_target_distribution) ); echo -n $(csg_get_interaction_property name),${tgt_dist}:')"
  upd_pots_flag="--upd-pots $do_potential_list"
  tgt_dists_flag="--tgt-dists $tgt_dist_list"
fi

# weather to make dU = 0 one point before the cut off
[[ "$(csg_get_property cg.inverse.iie.flatten_at_cut_off)" == 'true' ]] && flatten_at_cut_off_flag="--flatten-at-cut-off"

# RDF calculation
# for_all not necessary for most sim_prog, but also doesn't hurt.
for_all "non-bonded bonded" do_external rdf "$sim_prog"
# calculate distributions intramolecular
if [[ $(csg_get_property cg.inverse.iie.tgt_dcdh) != 'true' ]]; then
  for_all "non-bonded" do_external rdf "$sim_prog" --only-intra-nb
fi

# RDF extrapolation
g_extrap_factor=$(csg_get_property --allow-empty cg.inverse.iie.g_extrap_factor)
[[ -n $g_extrap_factor ]] && msg --color blue "Deprecated option g_extrap_factor will be ignored!"

# resample target distributions
for_all "non-bonded" do_external resample target --clean '$(csg_get_interaction_property inverse.target)' '$(csg_get_interaction_property name).dist.tgt'
# resample intramolecular only if needed and present. Later iie.py will only load the ones that are needed
if [[ $(csg_get_property cg.inverse.iie.tgt_dcdh) != 'true' ]]; then
  for_all "non-bonded" do_external resample target --no-extrap --skip-if-missing '$(csg_get_interaction_property inverse.target_intra)' '$(csg_get_interaction_property name).dist-intra.tgt'
fi

# Some arguments (cut_off, kBT) will be read directly from the settings.xml. They do not have a default in csg_defaults.xml.
# Others (closure, ...) could also be read from the settings file, but this bash script handles the defaults.
do_external update iie_pot "$iie_method" \
  ${verbose_flag-} \
  --closure "$(csg_get_property cg.inverse.iie.closure)" \
  --volume "$volume" \
  --topol "$topol" \
  --options "$CSGXMLFILE" \
  --g-tgt-ext "dist.tgt" \
  --g-cur-ext "dist.new" \
  ${g_intra_flag-} \
  --out "dpot.pure_iie" \
  ${pressure_constraint_flag-} \
  ${kirkwood_buff_constraint_flag-} \
  ${residual_weighting_flag-} \
  ${tgt_dcdh_flag-} \
  ${upd_pots_flag-} \
  ${tgt_dists_flag-} \
  ${flatten_at_cut_off_flag-}

# resample potentials. This is needed because non-bonded.max is sometimes larger than iie.cut-off and the former should define the end of the table
for_all "non-bonded" 'csg_resample --in $(csg_get_interaction_property name).dpot.pure_iie --out $(csg_get_interaction_property name).dpot.grid_adapted --grid $(csg_get_interaction_property min):$(csg_get_interaction_property step):$(csg_get_interaction_property max) --comment "adapted to grid in update_iie.sh"'
# csg_resample alone will sometimes lead to non-zero values, table extrapolate will make it zero
for_all "non-bonded" 'do_external table extrapolate --function constant --region right --no-flagupdate $(csg_get_interaction_property name).dpot.grid_adapted $(csg_get_interaction_property name).dpot.new'

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

for_all "bonded" do_external update ibi_single
