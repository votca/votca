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
${0##*/}, version 1.0_rc1 hgid: 49f54a9b1845112a273f8c1bf2c683f2674f71c7
This script:
- calculates potential for first parameter set (pot [function])
- flags that parameter set as 'active' (par pot)

Usage: ${0##*/}

USES: for_all csg_get_interaction_property do_external cp_from_main_dir

NEEDS: name
EOF
  exit 0
fi

check_deps "$0"

name=$(csg_get_interaction_property name)
property=$(csg_get_property cg.inverse.simplex.property)
function=$(csg_get_interaction_property inverse.simplex.function)
param_N=$(do_external pot $function --nparams)

# Calculate potential
do_external par pot simplex_$name.cur simplex_$name.new $param_N 0
