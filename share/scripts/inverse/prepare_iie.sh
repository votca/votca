#! /bin/bash
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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
This script prepares potentials in a generic way

Usage: ${0##*/}
EOF
   exit 0
fi

# case insensitive comparisons
shopt -s nocasematch

sim_prog="$(csg_get_property cg.inverse.program)"
method="$(csg_get_property cg.inverse.method)"
initial_guess_method="$(csg_get_property cg.inverse.${method}.initial_guess.method)"
main_dir=$(get_main_dir)
nb_interactions=$(csg_get_property --allow-empty cg.non-bonded.name)
kBT="$(csg_get_property cg.inverse.kBT)"
densities="$(csg_get_property cg.inverse.iie.densities)"
n_intra="$(csg_get_property cg.inverse.iie.n_intra)"
verbose=$(csg_get_property cg.inverse.iie.verbose)
cut_off="$(csg_get_property cg.inverse.iie.cut_off)"

if [ "${verbose}" == 'true' ]; then
    verbose_flag="--verbose"
elif [ "${verbose}" == 'false' ]; then
    verbose_flag=""
else
    die "verbose has to be 'true' or 'false'"
fi

case "$initial_guess_method" in
"table")
    for_all "bonded non-bonded" do_external prepare_single generic --use-table
    ;;
"bi")
    for_all "bonded non-bonded" do_external prepare_single generic --use-bi
    ;;
"ie")
    # only for IE
    initial_guess_closure="$(csg_get_property cg.inverse.${method}.initial_guess.closure)"
    initial_guess_ignore_intra="$(csg_get_property cg.inverse.${method}.initial_guess.ignore_intramolecular_correlation)"

    # bonded distributions inversion
    for_all "bonded" do_external prepare_single generic --use-bi
    # resample all target distributions
    for_all "non-bonded" do_external resample target '$(csg_get_interaction_property inverse.target)' '$(csg_get_interaction_property name).dist.tgt'

    # initial guess from rdf with hnc or py
    # TODO implement propper extrapolation (use raw, then extrapolate)
    #raw=""
    out=""
    for name in $nb_interactions; do
        #raw="$raw $(critical mktemp ${name}.pot.ie.raw.XXX)"
        out="$out ${name}.pot.new"
    done
    if [[ "${initial_guess_ignore_intra}" == "false" && $n_intra -gt 1 ]]; then
        critical cp -t . ${main_dir}/$(printf '%s.dist-incl.tgt' $nb_interactions)
        G_tgt_flag="--G-tgt $(printf '%s.dist-incl.tgt' $nb_interactions)"
    else
        G_tgt_flag=''
    fi

    closure_flag="--closure ${initial_guess_closure}"

    do_external dist invert_iie potential_guess \
    $verbose_flag \
    $closure_flag \
    --g-tgt $(printf "%s.dist.tgt" $nb_interactions) \
    $G_tgt_flag \
    --U-out $out \
    --kBT $kBT --densities $densities --cut-off $cut_off \
    --n-intra $n_intra
    ;;
*)
    die "cg.inverse.${method}.initial_guess has to be either table, bi, or ie"
    ;;
esac


if [[ $sim_prog != gromacs ]] ; then
  msg --color blue "######################################################"
  msg --color blue "# WARNING using this simulator is still experimental #"
  msg --color blue "# If you find a problem report it under:             #"
  msg --color blue "# https://github.com/votca/csg                       #"
  msg --color blue "######################################################"
fi
