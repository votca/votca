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

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
This script initizalizes potentials for simplex method 

Usage: ${0##*/}
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
#list of all parameters
parameters=( $(csg_get_property cg.non-bonded.inverse.simplex.parameters) )

for_all non-bonded do_external prepare_single simplex "${#parameters[@]}"

do_external simplex prepare_table "simplex.table.cur"
do_external simplex table_to_potentials "simplex.table.cur" "simplex.table.new"
echo "Transformation=None" > "simplex.state.new"

# cp confout.gro and so on
do_external prepare_generic $sim_prog

