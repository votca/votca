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
This script implements the update function for the Simplex Method.

Usage: ${0##*/}

USES:  die msg csg_get_interaction_property for_all do_external run_or_exit

NEEDS: cg.inverse.program
EOF
   exit 0
fi

check_deps "$0"

name=$(for_all non-bonded csg_get_interaction_property name);
param_N=$(do_external pot $function --nparams);

run_or_exit do_external update simplex_single

p_nr=$(grep -c 'pending$' simplex_$name.tmp);

if [ $p_nr == "0" ]; then
   # Generate new parameter set
   msg "Preparing new parameters"
   run_or_exit do_external update simplex_step simplex_$name.tmp simplex_$name.new $param_N
else 
   msg "Found 'pending' parameter set"
   cp simplex_$name.tmp simplex_$name.new
   cp state_$name.cur state_$name.new
fi