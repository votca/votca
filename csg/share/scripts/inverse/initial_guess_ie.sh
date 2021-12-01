#! /bin/bash
#
# Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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
This script calculates an integral equation initial guess

Usage: ${0##*/} [--help]
EOF
   exit 0
fi

verbose=$(csg_get_property cg.inverse.initial_guess.ie.verbose)
if [ "${verbose}" == 'true' ]; then
    verbose_flag="--verbose"
fi

main_dir=$(get_main_dir)
sim_prog="$(csg_get_property cg.inverse.program)"
nb_names=( $(for_all "non-bonded" csg_get_interaction_property name) )
nb_names="${nb_names[@]}"
ie_closure="$(csg_get_property cg.inverse.initial_guess.ie.closure)"
subtract_coulomb="$(csg_get_property cg.inverse.initial_guess.ie.subtract_coulomb)"
multistate="$(csg_get_property cg.inverse.multistate.enabled)"

# resample all target distributions
for_all "non-bonded" do_external resample target '$(csg_get_interaction_property inverse.target)' '$(csg_get_interaction_property name).dist.tgt'
for_all "non-bonded" do_external resample target --no-extrap '$(csg_get_interaction_property inverse.target_intra)' '$(csg_get_interaction_property name).dist-intra.tgt'

# topology for molecular conections and volume
topol=$(csg_get_property --allow-empty cg.inverse.initial_guess.ie.topol)
[[ -z $topol ]] && topol=$(csg_get_property cg.inverse.$sim_prog.topol)
[[ -f $topol ]] || die "${0##*/}: topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist"

# volume
volume=$(critical csg_dump --top "$topol" | grep 'Volume' | awk '{print $2}')
([[ -n "$volume" ]] && is_num "$volume") || die "could not determine the volume from file ${topol}"

# kBT
if [[ $multistate == true ]]; then
  state_nr=get_state_nr
  kbt=( $(csg_get_property cg.inverse.multistate.state_kBTs) )
  kbt="${kbt[$state_nr]}"
else
  kbt="$(csg_get_property cg.inverse.kBT)"
fi

# subtract Coulomb term
if [[ $subtract_coulomb == true ]]; then
  subtract_coulomb_flag="--subtract_coulomb"
fi


msg "Using initial guess for non-bonded interactions using integral equations"
# Some arguments (cut_off, cut_residual) will be read directly from the settings.xml. They do not have a default in csg_defaults.xml.
# Others (closure, ...) could also be read from the settings file, but this bash script handles the defaults.
do_external dist invert_iie potential_guess \
    ${verbose_flag-} \
    --closure "$ie_closure" \
    --volume "$volume" \
    --kBT "$kbt" \
    ${subtract_coulomb_flag-} \
    --topol "$topol" \
    --options "$CSGXMLFILE" \
    --g-tgt-ext "dist.tgt" \
    --g-tgt-intra-ext "dist-intra.tgt" \
    --out "pot.new"

# scale new potentials
scaling_factor_non_bonded="$(csg_get_property "cg.inverse.initial_guess.scale_non_bonded")"
if $(awk "BEGIN {exit (${scaling_factor_non_bonded} != 1.0 ? 0 : 1)}"); then
  for_all "non-bonded" 'mv $(csg_get_interaction_property name).pot.new $(csg_get_interaction_property name).pot.new.raw'
  for_all "non-bonded" 'do_external table linearop $(csg_get_interaction_property name).pot.new.raw $(csg_get_interaction_property name).pot.new '"${scaling_factor_non_bonded} 0"
fi

# resample potentials. This is needed because non-bonded.max is usually larger than iie.cut-off and the former should define the table
for_all "non-bonded" 'csg_resample --in $(csg_get_interaction_property name).pot.new --out $(csg_get_interaction_property name).pot.new --grid $(csg_get_interaction_property min):$(csg_get_interaction_property step):$(csg_get_interaction_property max) --comment "adapted to grid in initial_guess_ie.sh"'
# csg_resample alone will not do the job
for_all "non-bonded" 'do_external table extrapolate --function constant --region right $(csg_get_interaction_property name).pot.new $(csg_get_interaction_property name).pot.new'

# overwrite with .pot.in if table_overwrite and present
for_all "non-bonded" do_external prepare_single generic --table-overwrite
