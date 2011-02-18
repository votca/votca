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
This script implemtents the function initialize

Usage: ${0##*/}

Used external packages: gromacs
EOF
  exit 0
fi

from="$(csg_get_property cg.inverse.initial_configuration "laststep")"
conf="$(csg_get_property cg.inverse.gromacs.conf "conf.gro")"
if [ "$from" = "laststep" ]; then
  confout="$(csg_get_property cg.inverse.gromacs.conf_out "confout.gro")"
  cp_from_last_step "${confout}"
  critical mv "${confout}" "$conf"
elif [ "$from" = "maindir" ]; then
  cp_from_main_dir "$conf"
else
  die "${0##*/}: initial_configuration '$from' not implemented"
fi

mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
[ -f "$mdp" ] || die "${0##*/}: gromacs mdp file '$mdp' not found"

#convert potential in format for sim_prog
for_all non-bonded do_external convert_potential gromacs

check_temp "$mdp"
for_all "non-bonded" check_cutoff "$mdp"

grompp="$(csg_get_property cg.inverse.gromacs.grompp.bin "grompp")"
[ -n "$(type -p $grompp)" ] || die "${0##*/}: grompp binary '$grompp' not found"

index="$(csg_get_property cg.inverse.gromacs.grompp.index "index.ndx")"
[ -f "$index" ] || die "${0##*/}: grompp index file '$index' not found"
top="$(csg_get_property cg.inverse.gromacs.grompp.topol "topol.top")"
[ -f "$top" ] || die "${0##*/}: grompp topol file '$top' not found"
tpr="$(csg_get_property cg.inverse.gromacs.topol "topol.tpr")"
opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"

critical $grompp -n "${index}" -f "${mdp}" -p "$top" -o "$tpr" -c "${conf}" ${opts}
[ -f "$tpr" ] || die "${0##*/}: gromacs tpr file '$tpr' not found after runing grompp"
