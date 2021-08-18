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

Usage: ${0##*/}
EOF
   exit 0
fi

verbose=$(csg_get_property cg.inverse.initial_guess.ie.verbose)
if [ "${verbose}" == 'true' ]; then
    verbose_flag="--verbose"
elif [ "${verbose}" == 'false' ]; then
    verbose_flag=""
else
    die "verbose has to be 'true' or 'false'"
fi

main_dir=$(get_main_dir)
sim_prog="$(csg_get_property cg.inverse.program)"
nb_names=( $(for_all "non-bonded" csg_get_interaction_property name) )
nb_names="${nb_names[@]}"
kBT="$(csg_get_property cg.inverse.kBT)"
cut_off="$(csg_get_property cg.inverse.initial_guess.ie.cut_off)"
g_min="$(csg_get_property cg.inverse.initial_guess.ie.g_min)"
ie_closure="$(csg_get_property cg.inverse.initial_guess.ie.closure)"

# resample all target distributions
# TODO: resample dist-incl.tgt or better dist-intra.tgt
# TODO: one might want longer tgt RDF for initial guess but short for iterative update with extrapolation or vice-versa
for_all "non-bonded" do_external resample target '$(csg_get_interaction_property inverse.target)' '$(csg_get_interaction_property name).dist.tgt'

# initial guess from rdf with hnc or py
# TODO: implement propper extrapolation (use raw, then extrapolate)
for nb_name in $nb_names; do
  critical cp -t . "${main_dir}/${nb_name}.dist-incl.tgt"
done

# topology for molecular conections and volume
topol=$(csg_get_property --allow-empty cg.inverse.initial_guess.ie.topol)
[[ -z $topol ]] && topol=$(csg_get_property cg.inverse.$sim_prog.topol)
[[ -f $topol ]] || die "${0##*/}: topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist"

# volume
volume=$(critical csg_dump --top "$topol" | grep 'Volume' | awk '{print $2}')
([[ -n "$volume" ]] && is_num "$volume") || die "could not determine the volume from file ${topol}"

# do not put quotes around arguments with values ($G_tgt_flag)!
# this will give a codacy warning :/
msg "Using initial guess for non-bonded interactions using integral equations"
# Some arguments (cut_off, kBT) will be read directly from the settings.xml. They do not have a default in csg_defaults.xml.
# Others could also be read from the settings file, but this bash script handles the defaults.
do_external dist invert_iie potential_guess \
    "$verbose_flag" \
    --closure "$ie_closure" \
    --g-min "$g_min" \
    --volume "$volume" \
    --topol "$topol" \
    --options "$CSGXMLFILE" \
    --g-tgt-ext ".dist.tgt" \
    --G-tgt-ext ".dist-incl.tgt" \
    --U-out-ext ".pot.new"

# scale new potentials
scaling_factor_non_bonded="$(csg_get_property "cg.inverse.initial_guess.scale_non_bonded")"
if $(awk "BEGIN {exit (${scaling_factor_non_bonded} != 1.0 ? 0 : 1)}"); then
  for_all "non-bonded" 'mv $(csg_get_interaction_property name).pot.new $(csg_get_interaction_property name).pot.new.raw'
  for_all "non-bonded" 'do_external table linearop $(csg_get_interaction_property name).pot.new.raw $(csg_get_interaction_property name).pot.new '"${scaling_factor_non_bonded} 0"
fi

# overwrite with .pot.in if table_overwrite and present
for_all "non-bonded" do_external prepare_single generic --table-overwrite
