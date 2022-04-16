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

iie_algorithm="$(csg_get_property cg.inverse.iie.algorithm)"
sim_prog="$(csg_get_property cg.inverse.program)"

# target dc/dh
if [[ $(csg_get_property cg.inverse.iie.use_target_dcdh) == 'true' ]]; then
  tgt_dcdh_flag="--tgt-dcdh $(get_main_dir)/dcdh.npz"
else
  g_intra_flag="--g-cur-intra-ext dist-intra.new"
fi

# topology for molecular conections and volume
topol=$(csg_get_property cg.inverse.topol_xml)
[[ -f $topol ]] || die "${0##*/}: topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist"

# volume and k_B*T
volume=$(critical csg_dump --top "$topol" | grep 'Volume' | awk '{print $2}')
([[ -n "$volume" ]] && is_num "$volume") || die "could not determine the volume from file ${topol}."
kbt="$(csg_get_property cg.inverse.kBT)"
([[ -n "$kbt" ]] && is_num "$kbt") || die "could not determine kBT from options file."

# verbose
verbose=$(csg_get_property cg.inverse.verbose)
step_nr=$(get_current_step_nr)
[[ "${verbose}" == 'true' ]] && verbose_flag="--verbose"
[[ "${verbose}" == 'step0+1' ]] && [[ $step_nr == '0' || $step_nr == '1' ]] && verbose_flag="--verbose"

# RDF calculation
# for_all not necessary for most sim_prog, but also doesn't hurt.
for_all "non-bonded bonded" do_external rdf "$sim_prog"
# calculate distributions intramolecular
if [[ $(csg_get_property cg.inverse.iie.use_target_dcdh) != 'true' ]]; then
  for_all "non-bonded" do_external rdf "$sim_prog" --only-intra-nb
fi

# cut residual
if [[ $iie_algorithm == 'newton' ]]; then
  cut_residual="$(csg_get_property cg.inverse.newton.cut_off)"
  ([[ -n "$cut_residual" ]] && is_num "$cut_residual") || die "could not get cut-off (./inverse/newton/cut_off) from options file."
elif [[ $iie_algorithm == 'gauss-newton' ]]; then
  cut_residual="$(csg_get_property cg.inverse.gauss_newton.cut_residual)"
  ([[ -n "$cut_residual" ]] && is_num "$cut_residual") || die "could not get cut-residual (./inverse/gauss_newton/cut_residual) from options file."
else
  die "${0##*/}: value of cg.inverse.iie.algorithm must be newton or gauss-newton"
fi

# RDF extrapolation
g_extrap_factor=$(csg_get_property --allow-empty cg.inverse.iie.g_extrap_factor)
[[ -n $g_extrap_factor ]] && msg --color blue "Deprecated option g_extrap_factor will be ignored!"

# resample target distributions
for_all "non-bonded" do_external resample target --clean '$(csg_get_interaction_property inverse.target)' '$(csg_get_interaction_property name).dist.tgt'
# resample intramolecular only if needed and present. Later iie.py will only load the ones that are needed
if [[ $(csg_get_property cg.inverse.iie.use_target_dcdh) != 'true' ]]; then
  for_all "non-bonded" resample_intra_if_present
fi

# improve Jacobian in RDF onset region
if [[ $(csg_get_property cg.inverse.iie.improve_jacobian_onset) == "true" ]]; then
  improve_jacobian_onset_flag="--improve-jacobian-onset"
  onset_thresholds_flag="--onset-thresholds $(csg_get_property cg.inverse.iie.onset_thresholds)"
fi

# Some arguments (cut_off, ...) will be read directly from the settings.xml. They do not have a default in csg_defaults.xml.
# Others (closure, ...) could also be read from the settings file, but this bash script handles the defaults.
msg "Generating Jacobian with integral equation theory"
do_external generate iie_jacobian jacobian \
  ${verbose_flag-} \
  --closure "$(csg_get_property cg.inverse.iie.closure)" \
  --volume "$volume" \
  --kBT "$kbt" \
  --topol "$topol" \
  --options "$CSGXMLFILE" \
  --g-tgt-ext "dist.tgt" \
  --g-cur-ext "dist.new" \
  --cut-residual "$cut_residual" \
  ${tgt_dcdh_flag-} \
  ${g_intra_flag-} \
  ${improve_jacobian_onset_flag-} \
  ${onset_thresholds_flag-} \
  --out "jacobian.npz"

if [[ $iie_algorithm == 'newton' ]]; then
  # update nb with Newton
  do_external update newton
elif [[ $iie_algorithm == 'gauss-newton' ]]; then
  # update nb with Gauss-Newton
  do_external update gauss_newton
fi

for_all "bonded" do_external update ibi_single
