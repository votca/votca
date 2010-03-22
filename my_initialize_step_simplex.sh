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
This script implemtents the function 
initialize_step for the Simplex Method

Usage: ${0##*/}

USES:  get_step_nr get_main_dir for_all csg_get_interaction_property msg cp_from_last_step mv

NEEDS: -
EOF
   exit 0
fi

check_deps "$0"

main_dir=$(get_main_dir);
step_nr=$(get_current_step_nr);
init_steps=$(wc -l $main_dir/simplex.in | awk '{print ($1)}');

if [ $step_nr -le $[$init_steps] ]; then
   msg "Initialization step "$step_nr
else
   msg "Simplex step "$[$step_nr-$init_steps]
   for_all non-bonded 'cp_from_last_step state.new'
   for_all non-bonded 'mv state.new state.cur'
fi

# Get new pot from last step and make it current
for_all non-bonded 'cp_from_last_step $(csg_get_interaction_property name).pot.new'
for_all non-bonded 'mv $(csg_get_interaction_property name).pot.new $(csg_get_interaction_property name).pot.cur'

# Get new simplex table from last step and make it current
for_all non-bonded 'cp_from_last_step simplex.new'
for_all non-bonded 'mv simplex.new simplex.cur'