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
This script calcs the pressure for gromacs and writes it to outfile

Usage: ${0##*/} outfile

Used external packages: gromacs
EOF
   exit 0
fi

[[ -z $1 ]] && die "${0##*/}: Missing argument"

tpr="$(csg_get_property cg.inverse.gromacs.g_energy.topol "topol.tpr")"
[ -f "$tpr" ] || die "${0##*/}: Gromacs tpr file '$tpr' not found"

g_energy="$(csg_get_property cg.inverse.gromacs.g_energy.bin "g_energy")"
[ -n "$(type -p ${g_energy})" ] || die "${0##*/}: g_energy binary '$g_energy' not found"


opts="$(csg_get_property --allow-empty cg.inverse.gromacs.g_energy.opts)"

nsteps=$(get_simulation_setting nsteps)
dt=$(get_simulation_setting dt)
equi_time="$(csg_get_property cg.inverse.gromacs.equi_time 0)"
first_frame="$(csg_get_property cg.inverse.gromacs.first_frame 0)"

begin="$(awk -v dt=$dt -v frames=$first_frame -v eqtime=$equi_time 'BEGIN{print (eqtime > dt*frames ? eqtime : dt*frames) }')"

echo "Running ${g_energy}"
output=$(echo Pressure | critical ${g_energy} -b "${begin}" -s "${tpr}" ${opts})
echo "$output"
#the number pattern '-\?[0-9][^[:space:]]*[0-9]' is ugly, but it supports X X.X X.Xe+X Xe-X and so on
p_now=$(echo "$output" | sed -n 's/^Pressure[^-0-9]*\(-\?[0-9][^[:space:]]*[0-9]\)[[:space:]].*$/\1/p' ) || \
  die "${0##*/}: awk failed"
[ -z "$p_now" ] && die "${0##*/}: Could not get pressure from simulation"
echo "Pressure=${p_now}" > "$1"
