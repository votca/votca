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
This script implements the function update for the Simplex Method for a single interaction pair.

Usage: ${0##*/}

USES:  do_external die csg_get_property csg_get_interaction_property msg run_or_exit grep sed

NEEDS: name step min max inverse.do_potential
EOF
   exit 0
fi

check_deps "$0"

sim_prog=$(csg_get_property cg.inverse.program)
property=$(csg_get_property cg.inverse.simplex.property)
name=$(csg_get_interaction_property name);
function=$(csg_get_interaction_property inverse.simplex.function);
param_N=$(do_external pot $function --nparams | tail -1);


msg "Calc $property"
run_or_exit do_external $property $sim_prog

# For active parameter set, calculate ftar
if [ $(grep -c 'active$' simplex_$name.cur) == "1" ]; then
a_line_nr=$(($(grep -n -m1 'active$' simplex_$name.cur | sed 's/:.*//')-1));
else 
  die "Error: No 'active' parameter set found."
fi

msg "Calc $property ftar"
do_external update ftar_$property simplex_$name.cur simplex_$name.tmp $param_N $(($a_line_nr))
