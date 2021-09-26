#! /bin/bash
#
# Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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
This script calcs the pressure for lammps and writes it to outfile

Usage: ${0##*/} outfile

Used external packages: lammps
EOF
   exit 0
fi

[[ -z $1 ]] && die "${0##*/}: Missing argument"

p_file="$(csg_get_property cg.inverse.lammps.pressure_file)" 

[[ -f ${p_file} ]] || die "${0##*/}: pressure file '${p_file}' doesn't exist" 

p_now=$(awk 'NR > 1 {avg += $1} END {printf "%.16f\n", avg/(NR-1)}' ${p_file}) || die "${0##*/}: pressure averaging failed" 
echo "Pressure=${p_now}" > "$1"
