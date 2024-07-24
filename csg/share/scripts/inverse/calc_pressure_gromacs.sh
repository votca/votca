#! /usr/bin/env bash
#
# Copyright 2009-2015 The VOTCA Development Team (http://www.votca.org)
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

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
This script calcs the pressure for gromacs and writes it to outfile

Usage: ${0##*/} outfile

Used external packages: gromacs
EOF
   exit 0
fi

[[ -z $1 ]] && die "${0##*/}: Missing argument"

topol="$(csg_get_property --allow-empty cg.inverse.gromacs.g_energy.topol)"
[[ -z $topol ]] && topol=$(csg_get_property cg.inverse.gromacs.topol)
[[ -f $topol ]] || die "${0##*/}: Gromacs tpr file '$topol' not found"

g_energy=( $(csg_get_property cg.inverse.gromacs.g_energy.bin) )
[[ -n "$(type -p ${g_energy[0]})" ]] || die "${0##*/}: g_energy binary '${g_energy[0]}' not found"


opts="$(csg_get_property --allow-empty cg.inverse.gromacs.g_energy.opts)"

begin="$(calc_begin_time)"
if [[ ${CSG_RUNTEST} ]] && csg_calc "$begin" ">" "0"; then
  msg --color blue --to-stderr "Automatically setting begin time to 0, because CSG_RUNTEST was set"
  begin=0
fi

echo "Running ${g_energy[@]}"
#no critical here so that we can print the error
output=$(echo Pressure | ${g_energy[@]} -b "${begin}" -s "${topol}" ${opts} 2>&1)
ret="$?"
echo "$output" | gromacs_log "${g_energy[@]} -b "${begin}" -s "${topol}" ${opts}"
[[ $ret -eq 0 ]] || die "${0##*/}: '${g_energy[@]} -b "${begin}" -s "${topol}" ${opts}' failed"
#the number pattern '-\?[0-9][^[:space:]]*[0-9]' is ugly, but it supports X X.X X.Xe+X Xe-X and so on
#awk 'print $2' does not work for older version of g_energy as the format varies between
#^Pressure XXX (bar) and ^Pressure (bar) XXX
p_now=$(echo "$output" | sed -n 's/^Pressure[^-0-9]*\(\(-\?[0-9][^[:space:]]*[0-9]\|nan\)\)[[:space:]].*$/\1/p' ) || \
  die "${0##*/}: sed grep of Pressure failed"
[[ -z $p_now ]] && die "${0##*/}: Could not get pressure from simulation"

[[ $p_now = nan && $(csg_get_property cg.inverse.gromacs.g_energy.pressure.allow_nan) = no ]] && \
  die "${0##*/}: Pressure was nan, check your simulation (this usually means system has blow up -> use pre simulation)"

echo "Pressure=${p_now}" > "$1"
