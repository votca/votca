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
This script implements the potential update for the multistate iterative integral equation
methods

Usage: ${0##*/}  [--help]
EOF
   exit 0
fi

iie_algorithm="$(csg_get_property cg.inverse.iie.algorithm)"
sim_prog="$(csg_get_property cg.inverse.program)"
state_names="$(csg_get_property cg.inverse.multistate.state_names)"
state_names_arr=( $(csg_get_property cg.inverse.multistate.state_names) )
state_kBTs_arr=( $(csg_get_property cg.inverse.multistate.state_kBTs) )

# newton not supported any more
if [[ $iie_algorithm != gauss-newton ]]; then
  die "Multistate only works with the Gauss-Newton algorithm (more residuals than potentials)"
fi

# verbose
verbose=$(csg_get_property cg.inverse.verbose)
step_nr=$(get_current_step_nr)
[[ "${verbose}" == 'true' ]] && verbose_flag="--verbose"
[[ "${verbose}" == 'step0+1' ]] && [[ $step_nr == '0' || $step_nr == '1' ]] && verbose_flag="--verbose"

# cut residual
cut_residual="$(csg_get_property cg.inverse.gauss_newton.cut_residual)"
([[ -n "$cut_residual" ]] && is_num "$cut_residual") || die "could not get cut-residual (./inverse/gauss_newton/cut_residual) from options file."

# improve Jacobian in RDF onset region
if [[ $(csg_get_property cg.inverse.iie.improve_jacobian_onset) == "true" ]]; then
  improve_jacobian_onset_flag="--improve-jacobian-onset"
  onset_thresholds_flag="--onset-thresholds $(csg_get_property cg.inverse.iie.onset_thresholds)"
fi

# calc RDFs per state
# for_all not necessary for most sim_prog, but also doesn't hurt.
for state in $state_names; do
  if is_done "rdf_calculation_$state"; then
    continue
  fi
  pushd "$state"
  for_all "non-bonded bonded" do_external rdf $sim_prog
  # calculate distributions intramolecular
  if [[ $(csg_get_property cg.inverse.iie.use_target_dcdh) != 'true' ]]; then
    for_all "non-bonded" do_external rdf "$sim_prog" --only-intra-nb
  fi
  popd
  mark_done "rdf_calculation_$state"
done

# resample target distributions per state
for state in $state_names; do
  if is_done "resample_tgt_rdf_$state"; then
    continue
  fi
  pushd "$state"
  for_all "non-bonded" do_external resample target --clean '$(csg_get_interaction_property inverse.target)' '$(csg_get_interaction_property name).dist.tgt'
  if [[ $(csg_get_property cg.inverse.iie.use_target_dcdh) != 'true' ]]; then
    # resample intramolecular only if present. Later iie.py will only load the ones that are needed
    for_all "non-bonded" resample_intra_if_present
  fi
  popd
  mark_done "resample_tgt_rdf_$state"
done

# calculate Jacobian per state
for s in "${!state_names_arr[@]}"; do
  state="${state_names_arr[s]}"
  kBT="${state_kBTs_arr[s]}"
  if is_done "jacobian_$state";
  then
    continue
  fi
  pushd "$state"
  msg "calculating Jacobian for state $state"

  # target dc/dh
  if [[ $(csg_get_property cg.inverse.iie.use_target_dcdh) == true ]]; then
    tgt_dcdh_flag="--tgt-dcdh $(get_main_dir)/${state}/dcdh.npz"
  else
    g_intra_flag="--g-cur-intra-ext dist-intra.new"
  fi

  # topology and volume
  topol_state="$(csg_get_property cg.inverse.topol_xml)"
  [[ -f $topol_state ]] || die "${0##*/}: topol file '$topol_state' not found, possibly you have to add it to cg.inverse.filelist"
  volume_state=$(critical csg_dump --top "$topol_state" | grep 'Volume' | awk '{print $2}')
  ([[ -n "$volume_state" ]] && is_num "$volume_state") || die "could not determine the volume from file ${topol_state}"

  # Some arguments (cut_off, ...) will be read directly from the settings.xml. They do not have a default in csg_defaults.xml.
  # Others (closure, ...) could also be read from the settings file, but this bash script handles the defaults.
  do_external generate iie_jacobian jacobian \
    ${verbose_flag-} \
    --closure "$(csg_get_property cg.inverse.iie.closure)" \
    --volume $volume_state \
    --kBT "$kBT" \
    --topol $topol_state \
    --options "$CSGXMLFILE" \
    --g-tgt-ext "dist.tgt" \
    --g-cur-ext "dist.new" \
    --cut-residual "$cut_residual" \
    ${tgt_dcdh_flag-} \
    ${g_intra_flag-} \
    ${improve_jacobian_onset_flag-} \
    ${onset_thresholds_flag-} \
    --out "jacobian.npz"
  popd
  mark_done "jacobian_$state";
done

# update nb with Gauss-Newton
do_external update multistate_gauss_newton

# TODO bonded multistate average
for_all "bonded" do_external update ibi_single
