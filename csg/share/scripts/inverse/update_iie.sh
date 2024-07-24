#! /usr/bin/env bash
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
This script implements the function update for the HNC methods

There is not support for not doing potential updates on some interactions yet.
Also no support for per interaction extrapolation.


Usage: ${0##*/}
EOF
   exit 0
fi

iie_closure="$(csg_get_property cg.inverse.iie.closure)"
iie_method="$(csg_get_property cg.inverse.iie.method)"
ignore_intramolecular_correlation="$(csg_get_property cg.inverse.iie.ignore_intramolecular_correlation)"
sim_prog="$(csg_get_property cg.inverse.program)"
nb_interactions=$(csg_get_property --allow-empty cg.non-bonded.name)

# TODO: outsource the extrapolation!
if [[ $iie_method == 'gauss-newton' ]]; then
    pressure_constraint=$(csg_get_property cg.inverse.iie.pressure_constraint)
    if is_num "${pressure_constraint}"; then
        p_file="${name}.pressure"
        do_external pressure "$sim_prog" "$p_file"
        p_now="$(sed -n 's/^Pressure=\(.*\)/\1/p' "$p_file")" || die "${0##*/}: sed of Pressure failed"
        [[ -z $p_now ]] && die "${0##*/}: Could not get pressure from simulation"
        echo "New pressure $p_now, target pressure $pressure_constraint"
        pressure_constraint_flag="--pressure-constraint $pressure_constraint,$p_now"
    else
        pressure_constraint_flag=""
    fi
    extrap_near_core_flag="--extrap-near-core $(csg_get_property cg.inverse.iie.extrap_near_core)"
    fix_near_cut_off_flag="--fix-near-cut-off $(csg_get_property cg.inverse.iie.fix_near_cut_off)"
elif [[ $iie_method == 'newton' || $iie_method == 'newton-mod'  ]]; then
    if [[ $(csg_get_property cg.inverse.iie.cut_jacobian) == 'true' ]]; then
        cut_jacobian_flag="--cut-jacobian"
    fi
fi

kBT="$(csg_get_property cg.inverse.kBT)"
densities="$(csg_get_property cg.inverse.iie.densities)"
n_intra="$(csg_get_property cg.inverse.iie.n_intra)"
verbose=$(csg_get_property cg.inverse.iie.verbose)
cut_off="$(csg_get_property cg.inverse.iie.cut_off)"
g_extrap_factor="$(csg_get_property cg.inverse.iie.g_extrap_factor)"

if [[ "${verbose}" == 'true' ]]; then
    verbose_flag="--verbose"
elif [[ "${verbose}" == 'false' ]]; then
    verbose_flag=""
else
    die "verbose has to be 'true' or 'false'"
fi

if [[ "${ignore_intramolecular_correlation}" == 'false' ]]; then
    # TODO: n_intra will be an array, check if any greater one
    if [[ $n_intra -gt 1 ]]; then
        # calculate distributions including intra
        for_all "non-bonded" do_external rdf_incl_intra generic --include-intra
        G_cur_flag="--G-cur $(printf "%s.dist-incl.new" "$nb_interactions")"
    else
        G_cur_flag=""
    fi
elif [[ "${ignore_intramolecular_correlation}" == 'true' ]]; then
    G_cur_flag=""
else
    die "ignore_intramolecular_correlation has to be 'true' or 'false'"
fi

# for_all not necessary for most sim_prog, but also doesn't hurt.
for_all "non-bonded bonded" do_external rdf "$sim_prog"
for_all "non-bonded" do_external resample target '$(csg_get_interaction_property inverse.target)' '$(csg_get_interaction_property name).dist.tgt'

# do not put quotes around arguments with values ($G_cur_flag and some others)!
# this will give a codacy warning :/
do_external update iie_pot "$iie_method" \
    "$verbose_flag" \
    --closure "$iie_closure" \
    --g-tgt $(printf "%s.dist.tgt" "$nb_interactions") \
    --g-cur $(printf "%s.dist.new" "$nb_interactions") \
    $G_cur_flag \
    --U-cur $(printf "%s.pot.cur" "$nb_interactions") \
    --U-out $(printf "%s.dpot.new" "$nb_interactions") \
    --kBT "$kBT" --densities "$densities" --cut-off "$cut_off" \
    --g-extrap-factor "$g_extrap_factor" \
    $extrap_near_core_flag \
    $fix_near_cut_off_flag \
    $pressure_constraint_flag \
    $cut_jacobian_flag \
    --n-intra "$n_intra"

for_all "bonded" do_external update ibi_single
