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
This script runs gromacs
for the Inverse Boltzmann Method

Usage: ${0##*/}

USES: run_or_exit use_mpi csg_get_property check_deps

OPTIONAL: cg.inverse.mpi.cmd cg.inverse.gromacs.mdrun.opts cg.inverse.gromacs.topol cg.inverse.gromacs.traj_type cg.inverse.gromacs.mdrun.bin
EOF
   exit 0
fi

tpr="$(csg_get_property cg.inverse.gromacs.topol "topol.tpr")"
[ -f "$tpr" ] || die "${0##*/}: gromacs tpr file '$tpr' not found"

mdrun="$(csg_get_property cg.inverse.gromacs.mdrun.bin "mdrun")"
[ -n "$(type -p $mdrun)" ] || die "${0##*/}: mdrun binary '$mdrun' not found"

opts="$(csg_get_property --allow-empty cg.inverse.gromacs.mdrun.opts)"

check_deps "$0"

if use_mpi; then
  mpicmd=$(csg_get_property --allow-empty cg.inverse.mpi.cmd)
  run_or_exit $mpicmd $mdrun -s "${tpr}" ${opts}
else
  run_or_exit $mdrun -s "${tpr}" ${opts}
fi

ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
traj="traj.${ext}"
[ -f "$traj" ] || die "${0##*/}: gromacs traj file '$traj' not found after mdrun"
