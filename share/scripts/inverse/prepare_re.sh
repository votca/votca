#! /bin/bash
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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
This script implements the preparation of the relative entropy method iteration

Usage: ${0##*/}
EOF
   exit 0
fi


# get initial parameters from main dir and make it current parameters
for_all non-bonded 'cp_from_main_dir --rename $(csg_get_interaction_property name).param.init $(csg_get_interaction_property name).param.new'

sim_prog="$(csg_get_property cg.inverse.program)"
#cp confout.gro and so on
do_external prepare_generic $sim_prog

# run csg_reupdate to generate intital potential tables
msg --color green "Generating potential tables from the initial parameters"
critical csg_reupdate --gentable true --param-in-ext param.new --options $CSGXMLFILE 
