#! /usr/bin/env bash
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
method="$(csg_get_property cg.inverse.method)"

for_all "bonded non-bonded" do_external prepare_single $method

if [[ $sim_prog != @(gromacs|lammps) ]] ; then
  msg --color blue "######################################################"
  msg --color blue "# WARNING using this simulator is still experimental #"
  msg --color blue "# If you find a problem report it under:             #"
  msg --color blue "# https://github.com/votca/votca/issues              #"
  msg --color blue "######################################################"
fi

dihedral_names=$(for_all "dihedral" csg_get_interaction_property name)
if [[ ${dihedral_names} ]]; then
  msg --color blue "###########################################################"
  msg --color blue "# WARNING using VOTCA for dihedrals is still experimental #"
  msg --color blue "# If you find any issues please report them at:           #"
  msg --color blue "# https://github.com/votca/votca/issues                   #"
  msg --color blue "###########################################################"
fi
