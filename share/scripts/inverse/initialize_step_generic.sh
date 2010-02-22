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
for the Inverse Boltzmann Method

Usage: ${0##*/}

USES:  for_all csg_get_interaction_property mv check_deps cp_from_last_step

NEEDS: name 
EOF
   exit 0
fi

check_deps "$0"

#get new pot from last step 
for_all non-bonded 'cp_from_last_step $(csg_get_interaction_property name).pot.new'

#make it current potential 
for_all non-bonded 'mv $(csg_get_interaction_property name).pot.new $(csg_get_interaction_property name).pot.cur'

