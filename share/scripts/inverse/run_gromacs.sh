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
This script runs a gromacs simulation

Usage: ${0##*/}

Used external packages: gromacs
EOF
   exit 0
fi

tpr="$(csg_get_property cg.inverse.gromacs.topol "topol.tpr")"


mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
[ -f "$mdp" ] || die "${0##*/}: gromacs mdp file '$mdp' not found"

conf="$(csg_get_property cg.inverse.gromacs.conf "conf.gro")"
[ -f "$conf" ] || die "${0##*/}: gromacs initial configuration file '$conf' not found"

confout="$(csg_get_property cg.inverse.gromacs.conf_out "confout.gro")"
mdrun_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.mdrun.opts)"

index="$(csg_get_property cg.inverse.gromacs.grompp.index "index.ndx")"
[ -f "$index" ] || die "${0##*/}: grompp index file '$index' not found"
top="$(csg_get_property cg.inverse.gromacs.grompp.topol "topol.top")"
[ -f "$top" ] || die "${0##*/}: grompp topol file '$top' not found"

grompp_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"

grompp="$(csg_get_property cg.inverse.gromacs.grompp.bin "grompp")"
[ -n "$(type -p $grompp)" ] || die "${0##*/}: grompp binary '$grompp' not found"

critical $grompp -n "${index}" -f "${mdp}" -p "$top" -o "$tpr" -c "${conf}" ${grompp_opts}
[ -f "$tpr" ] || die "${0##*/}: gromacs tpr file '$tpr' not found after runing grompp"

mdrun="$(csg_get_property cg.inverse.gromacs.mdrun.bin "mdrun")"
#no check for mdrun, because mdrun_mpi could maybe exist only computenodes

mpicmd=$(csg_get_property --allow-empty cg.inverse.parallel.cmd)
critical $mpicmd $mdrun -s "${tpr}" -c "${confout}" -o traj.trr -x traj.xtc ${mdrun_opts}

ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
traj="traj.${ext}"
[ -f "$traj" ] || die "${0##*/}: gromacs traj file '$traj' not found after running mdrun"

[ -f "$confout" ] || die "${0##*/}: Gromacs end coordinate '$confout' not found after running mdrun"
