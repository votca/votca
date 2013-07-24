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
This script runs lammps for the Inverse Boltzmann Method

Usage: ${0##*/}

Used external packages: lammps
EOF
    exit 0
fi

script="$(csg_get_property cg.inverse.lammps.script)"
[[ -f $script ]] || die "${0##*/}: lammps script '$script' not found (make sure it is in cg.inverse.filelist)"

lmp_bin="$(csg_get_property cg.inverse.lammps.command)"
#no check for Espresso, because Espresso could maybe exist only computenodes

traj="$(csg_get_property cg.inverse.lammps.traj)"

if [[ -n $CSGENDING ]]; then
  echo "${0##*/} does not support wallclock time yet (go here and implement it). Per step wallclock time check is still performed!"
fi

method="$(csg_get_property cg.inverse.method)"
[[ $method = "ibi" ||  $method = imc ]] || die "${0##*/}: lammps only supports method: ibi"

critical $lmp_bin -in "${script}"

[[ -f $traj ]] || die "${0##*/}: traj file '$traj' wasn't found after running lammps"
