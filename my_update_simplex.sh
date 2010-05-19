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
- implements the update function for each non-bonded interaction
- performs downhill simplex algorithm if no pending parameter sets present
- continues with next parameter set in table if otherwise

Usage: ${0##*/}

USES:  for_all csg_get_interaction_property do_external run_or_exit

NEEDS: name function
EOF
   exit 0
fi

check_deps "$0"

name=$(for_all non-bonded csg_get_interaction_property name);
function=$(for_all non-bonded csg_get_interaction_property inverse.simplex.function);
param_N=$(do_external pot $function --nparams);

run_or_exit for_all non-bonded do_external update simplex_single

p_nr=$(grep -c 'pending$' simplex_$name.tmp);

if [ $p_nr == "0" ]; then
   # Generate new parameter set
   msg "Calculating new parameter set"
   run_or_exit do_external update simplex_step simplex_$name.tmp simplex_$name.new $param_N
else 
   msg "Continuing with next parameter set"
   cp simplex_$name.tmp simplex_$name.new
   cp state_$name.cur state_$name.new
fi