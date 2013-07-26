#! /bin/bash
#
# Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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
This script runs a generic simulation program

Usage: ${0##*/}
EOF
    exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
script="$(csg_get_property cg.inverse.$sim_prog.script)"
[[ -f $script ]] || die "${0##*/}: $sim_prog script '$script' not found (make sure it is in cg.inverse.filelist)"

cmd="$(csg_get_property cg.inverse.$sim_prog.command)"
#no check for Espresso, because Espresso could maybe exist only computenodes

opts=$(csg_get_property --allow-empty cg.inverse.$sim_prog.opts)
#expand ${script} in there
opts="$(eval echo $opts)"

# trajectory file
traj="$(csg_get_property cg.inverse.$sim_prog.traj)"

if [[ -n $CSGENDING ]]; then
  echo "${0##*/}: $sim_prog does not support wallclock time yet (go here and implement it). Per step wallclock time check is still performed!"
fi

method="$(csg_get_property cg.inverse.method)"
shopt -s extglob
[[ $method = @(ibi|imc|optimizer|re) ]] || die "${0##*/}: ${sim_prog} does not support method $method yet!"

critical $cmd ${opts}

simulation_finish
