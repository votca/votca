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
- finds the pending parameter set
- calculates the potential for this set (pot [function])
- flags this parameter set as 'active' (par pot)

Usage: ${0##*/} method

USES: for_all csg_get_interaction_property do_external run_or_exit

NEEDS: name function
EOF
   exit 0
fi

check_deps "$0"

name=$(for_all non-bonded csg_get_interaction_property name);
function=$(for_all non-bonded csg_get_interaction_property inverse.simplex.function);
param_N=$(do_external pot $function --nparams);
p_nr=$(grep -c 'pending$' simplex_$name.new);
tmp=$(mktemp simplex_XXX);

# Finds pending parameter set
if [ $p_nr > "0" ]; then
p_line_nr=$(($(grep -n -m1 'pending$' simplex_$name.new | sed 's/:.*//')-1));
else 
p_line_nr=0;
fi

# Calculate potential
for_all "non-bonded" \
  run_or_exit do_external par pot simplex_$name.new simplex_$name.new $param_N $p_line_nr