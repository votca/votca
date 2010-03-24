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
This script implemtents the function update
for the Simplex Method

Usage: ${0##*/}

USES:  die msg csg_get_property for_all do_external

NEEDS: cg.inverse.program
EOF
   exit 0
fi

check_deps "$0"

run_or_exit do_external update simplex_single

p_nr=$(grep -c 'pending$' simplex.tmp);

if [ $(grep -c 'pending$' simplex.cur) == "0" ]; then
   c_line_nr=$(($(grep -n -m1 'pending' simplex.tmp | sed 's/:.*//')-2));
   # Generate new parameter set
   msg "Preparing new parameters"
   run_or_exit do_external update simplex_step simplex.tmp simplex.new $p_nr
else 
   msg "Found 'pending' parameter set"
   cp state.cur state.new
fi

