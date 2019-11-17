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
This script initializes potentials for ihnc

Usage: ${0##*/}
EOF
   exit 0
fi

names=( $(csg_get_interaction_property --all name) )
if [[ ${#names[@]} -gt 1 ]]; then
    die "IHNC multicomponent is not yet implemented"
fi

[[ -n $(csg_get_property --allow-empty cg.bonded.name) ]] && die "IHNC does not support bonded interactions, go and implement it"

method="$(csg_get_property cg.inverse.method)"
for_all "non-bonded" do_external prepare_single $method
