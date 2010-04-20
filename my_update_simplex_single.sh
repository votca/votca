#! /bin/bash
# 
# Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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
for the Simplex Method for a single pair

Usage: ${0##*/} step_nr

USES:  do_external die csg_get_interaction_property msg run_or_exit wc awk

NEEDS: name step min max inverse.do_potential
EOF
   exit 0
fi

check_deps "$0"

sim_prog=$(csg_get_property cg.inverse.program)
property=$(csg_get_property cg.inverse.simplex.property)

#if [ $property eq "rdf" ]; then
msg "Calc rdf"
for_all non-bonded do_external rdf $sim_prog
#fi

if [ $property eq "surften"]; then
msg "Calc surften"
for_all non-bonded do_external surften $sim_prog
fi

# For active parameter set, calculate ftar
if [ $(grep -c 'active$' simplex.cur) == "1" ]; then
a_line_nr=$(($(grep -n -m1 'active$' simplex.cur | sed 's/:.*//')-1));
else 
  die "Error: No 'active' parameter set found."
fi

msg "Calc ftar"
for_all non-bonded do_external update simplex_ftar '$(csg_get_interaction_property name).dist.tgt' \
'$(csg_get_interaction_property name).dist.new' simplex.cur simplex.tmp $(($a_line_nr))
