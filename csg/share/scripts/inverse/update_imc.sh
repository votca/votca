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
This script implements the function update
for the Inverse Monte Carlo Method

Usage: ${0##*/}
EOF
   exit 0
fi

imc_algorithm="$(csg_get_property cg.inverse.imc.algorithm)"
sim_prog="$(csg_get_property cg.inverse.program)"

# IMC groups
imc_groups=$(csg_get_interaction_property --all inverse.imc.group)
imc_groups=$(remove_duplicate $imc_groups)
[[ -z ${imc_groups} ]] && die "${0##*/}: No imc groups defined"

# Regularization
default_reg=$(csg_get_property cg.inverse.imc.default_reg)
is_num "${default_reg}" || die "${0##*/}: value of cg.inverse.imc.default_reg should be a number"

# improve Jacobian in RDF onset region
if [[ $(csg_get_property cg.inverse.imc.improve_jacobian_onset) == "true" ]]; then
  improve_jacobian_onset_flag="--improve-jacobian-onset"
  onset_threshold_flag="--onset-threshold $(csg_get_property cg.inverse.imc.onset_threshold)"
fi

# old IMC Newton algorithm, using C++ for solving
# no constraints, but regularization is implemented
if [[ $imc_algorithm == 'newton' ]]; then

  [[ -n $improve_jacobian_onset_flag ]] && die "${0##*/}: for Newton IMC, the Jacobian improvement is not implemented."

  # calculate IMC matrix
  do_external imc_stat $sim_prog

  # IMC updates
  for group in $imc_groups; do
    if [[ $group == "none" ]]; then
      continue
    fi
    reg="$(csg_get_property cg.inverse.imc.${group}.reg ${default_reg})" #filter me away
    is_num "${reg}" || die "${0##*/}: value of cg.inverse.imc.${group}.reg should be a number"
    msg "solving linear equations for imc group '$group' (regularization ${reg})"
    critical csg_imc_solve --imcfile "${group}.imc" --gmcfile "${group}.gmc" --idxfile "${group}.idx" --regularization "${reg}"
  done

  for_all "non-bonded bonded" do_external update imc_single

# new IMC Gauss-Newton, using Python for solving
# constraints and weighting, but no regularization implemented
elif [[ $imc_algorithm == 'gauss-newton' ]]; then

  # no imc groups implemented for Gauss-Newton
  imc_groups=$(csg_get_interaction_property --all inverse.imc.group)
  imc_groups=$(remove_duplicate $imc_groups)
  imc_groups_array=( $imc_groups )
  [[ ${#imc_groups_array[@]} != 1 ]] && die "${0##*/}: for Gauss-Newton IMC, there can only be one IMC group"
  [[ ${imc_groups_array[0]} == none ]] && die "${0##*/}: only IMC group is none, needs to be something else"
  imc_group=${imc_groups[0]}

  # neither regularization
  [[ $default_reg != "0" ]] && die "${0##*/}: for Gauss-Newton IMC, regularization is not implemented."

  # topology for molecular conections and volume
  topol=$(csg_get_property cg.inverse.topol_xml)
  [[ -f $topol ]] || die "${0##*/}: topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist"

  # volume
  volume=$(critical csg_dump --top "$topol" | grep 'Volume' | awk '{print $2}')
  ([[ -n "$volume" ]] && is_num "$volume") || die "${0##*/}: could not determine the volume from file ${topol}"

  # verbose
  verbose=$(csg_get_property cg.inverse.verbose)
  step_nr=$(get_current_step_nr)
  [[ "${verbose}" == 'true' ]] && verbose_flag="--verbose"
  [[ "${verbose}" == 'step0+1' ]] && [[ $step_nr == '0' || $step_nr == '1' ]] && verbose_flag="--verbose"

  # check improve dist
  improve_dist_new_all=$(csg_get_interaction_property --all improve_dist_near_core.new)
  improve_dist_target_all=$(csg_get_interaction_property --all improve_dist_near_core.target)
  improve_dist_new_all=$(remove_duplicate $improve_dist_new_all)
  improve_dist_target_all=$(remove_duplicate $improve_dist_target_all)
  if [[ $(wc -w <<< "$improve_dist_new_all") > 1 || $(wc -w <<< "$improve_dist_target_all") > 1 ]]; then
    die "${0##*/}: If using IMC Gauss-Newton either all or no RDFs have to use improve_dist_near_core.{new,target}"
  fi

  # check that all nb interactions have same min/max/step
  check_same_grid() {
    for prop in min max step; do
      value=$(csg_get_interaction_property ${prop})
      prev_var="prev_${prop}"
      if [[ -n ${!prev_var} && ${!prev_var} != $value ]]; then
        die "${0##*/}: If using IMC Gauss-Newton all RDFs min, max, step have to be the same."
      fi
      declare prev_$prop=$value
    done
  }
  export -f check_same_grid
  for_all "non-bonded" check_same_grid

  # IMC matrix, RDF, target RDF
  use_target_matrix=$(csg_get_property cg.inverse.imc.use_target_matrix)
  if [[ "${use_target_matrix}" == 'true' ]]; then
    imc_group="$(get_main_dir)/$imc_group"
    [[ ( ! -f "${imc_group}.gmc" ) || ( ! -f "${imc_group}.idx" ) ]] && die "${0##*/}: if inverse.imc.use_target_matrix is true, .gmc and .idx file need to be in main dir."
    # calculate distributions only
    for_all "non-bonded bonded" do_external rdf "$sim_prog"
    # resample target distributions
    for_all "non-bonded" do_external resample target --clean '$(csg_get_interaction_property inverse.target)' '$(csg_get_interaction_property name).dist.tgt'
  else
    # calculate IMC matrix, RDF, and resample target
    do_external imc_stat $sim_prog
    if [[ $improve_dist_new_all == "true" ]]; then
      # raw current rdf needed for IMC onset fix
      g_raw_ext_flag="--g-cur-raw-ext=dist.new.raw"
    fi
  fi

  # convert IMC matrix to dg/du Jacobian
  msg "Converting imc matrix into dg/du Jacobian"
  do_external convert imc_matrix \
    ${verbose_flag-} \
    --options "$CSGXMLFILE" \
    --imc-matrix "${imc_group}.gmc" \
    --imc-index "${imc_group}.idx" \
    --volume "$volume" \
    --topol "${topol}" \
    --g-tgt-ext "dist.tgt" \
    --g-cur-ext "dist.new" \
    ${improve_jacobian_onset_flag-} \
    ${onset_threshold_flag-} \
    --out "jacobian.npz"

  # update nb with Gauss-Newton
  do_external update gauss_newton

  # do bonded (will in most cases do IBI)
  for_all "bonded" do_external update imc_single
else
  die "${0##*/}: value of cg.inverse.imc.algorithm must be newton or gauss-newton"
fi
