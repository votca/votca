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
This script:
- calculates the new property
- compares it to the target property and calculates the target function accordingly

Usage: ${0##*/}

USES:  csg_get_property csg_get_interaction_property do_external run_or_exit

NEEDS: name program property function
EOF
   exit 0
fi

check_deps "$0"

name=$(csg_get_interaction_property name)
sim_prog=$(csg_get_property cg.inverse.program)
property=$(csg_get_property cg.inverse.simplex.property)
count=$(echo "$property" | wc -w)
function=$(csg_get_interaction_property inverse.simplex.function)
param_N=$(do_external pot $function --nparams | tail -1)

# Calculate each property
for i in $property; do
  msg "Calc $i"
  do_external $i $sim_prog
  # Find 'active' parameter set
  if [ $(grep -c 'active$' simplex_$name.cur) == "1" ]; then
    a_line_nr=$(($(grep -n -m1 'active$' simplex_$name.cur | sed 's/:.*//')-1));
  else
    die "Error: No 'active' parameter set found."
  fi
done

# Calculate penalty function
if [ $count -eq "1" ]; then
 msg "Calc $property ftar"
 do_external update ftar_$property simplex_$name.cur simplex_$name.tmp $param_N $(($a_line_nr))
elif [ $count -gt "1" ]; then
  for i in $property; do
    msg "Calc $i ftar"
    do_external update ftar_$i simplex_$name.cur simplex_$name\_$i.tmp $param_N $(($a_line_nr))
  done
  msg "Calc total ftar"
  do_external update ftar_merge simplex_$name.tmp $param_N
fi
