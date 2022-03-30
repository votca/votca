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
This script prepares potentials in a generic way

Usage: ${0##*/}
EOF
  exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
initial_guess_method="$(csg_get_property cg.inverse.initial_guess.method)"
state_names="$(csg_get_property cg.inverse.multistate.state_names)"

case "$initial_guess_method" in
"table")
  for_all "bonded non-bonded" do_external prepare_single generic --use-table
  ;;
"bi"|"ie")
  # do potential guess for each state
  for state in $state_names; do
    pushd $state
    msg "for state ${state}:"
    # Boltzmann inverse
    if [[ $initial_guess_method == bi ]]; then
      for_all "bonded non-bonded" do_external prepare_single generic --use-bi
    # integral equation
    elif [[ $initial_guess_method == ie ]]; then
      bonded_interactions=( $(csg_get_property --allow-empty cg.bonded.name) )
      if [[ -n $bonded_interactions ]]; then
        for_all "bonded" do_external prepare_single generic --use-bi
      fi
      do_external initial_guess ie
    fi
    popd
  done
  # now average all states potential, not using weights, same as in MS-IBI Moore 2014
  export state_names  # needed for the following command
  for_all "non-bonded bonded" do_external table average --clean --output '$(csg_get_interaction_property name).pot.with_error' \
    '$(for s in $state_names; do echo $s/$(csg_get_interaction_property name).pot.new; done)'
  for_all "non-bonded bonded" do_external table remove_error '$(csg_get_interaction_property name).pot.with_error' '$(csg_get_interaction_property name).pot.new'
  ;;
*)
  die "cg.inverse.initial_guess.method has to be either table, bi, or ie"
  ;;
esac

if [[ $sim_prog != @(gromacs|lammps) ]] ; then
  msg --color blue "######################################################"
  msg --color blue "# WARNING using this simulator is still experimental #"
  msg --color blue "# If you find a problem report it under:             #"
  msg --color blue "# https://github.com/votca/csg                       #"
  msg --color blue "######################################################"
fi

dihedral_names=$(for_all "dihedral" csg_get_interaction_property name)
if [[ ${dihedral_names} ]]; then
  msg --color blue "###########################################################"
  msg --color blue "# WARNING using VOTCA for dihedrals is still experimental #"
  msg --color blue "# If you find any issues please report them at:           #"
  msg --color blue "# https://github.com/votca/csg                            #"
  msg --color blue "###########################################################"
fi