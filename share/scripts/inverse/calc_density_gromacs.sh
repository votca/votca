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

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
This script calcs the density for gromacs

Usage: ${0##*/} outputfile csg_density_options
EOF
   exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing argument"
output="$1"
shift

sim_prog="$(csg_get_property cg.inverse.program)"

if [[ $sim_prog = "gromacs" ]]; then
else
  die "${0##*/}: Simulation program '$sim_prog' not supported yet"
fi

topol=$(csg_get_property cg.inverse.gromacs.topol)
[[ -f $topol ]] || die "${0##*/}: gromacs topol file '$topol' not found"

traj=$(csg_get_property cg.inverse.$sim_prog.traj)
[[ -f $traj ]] || die "${0##*/}: traj file '$traj' not found"

name=$(csg_get_interaction_property name)

equi_time="$(csg_get_property cg.inverse.gromacs.equi_time)"
first_frame="$(csg_get_property cg.inverse.gromacs.first_frame)"

with_errors=$(csg_get_property cg.inverse.gromacs.density.with_errors)
if [[ ${with_errors} = "yes" ]]; then
  suffix="_with_errors"
else
  suffix=""
fi

if is_done "${name}_density_analysis${suffix}"; then
  echo "density analysis is already done"
  exit 0
fi

if [[ ${with_errors} = "yes" ]]; then
  msg "Calculating density for $name with errors"
  block_length=$(csg_get_property cg.inverse.gromacs.density.block_length)
  critical csg_density --trj "$traj" --top "$topol" --out "${output}.block" --begin "$equi_time" --first-frame "$first_frame" --block-length $block_length "$@"
  for i in ${output}.block_*; do
    [[ -f $i ]] || die "${0##*/}: Could not find ${output}.block_* after running csg_density, that usually means the blocksize (cg.inverse.gromacs.density.block_length) is too big."
  done
  #mind the --clean option to avoid ${name}.dist.block_* to fail on the second run
  do_external table average --clean --output "${output}" ${output}.block_*
else
  msg "Calculating density for $name"
  critical csg_density --trj "$traj" --top "$topol" --out "$output" --begin "$equi_time" --first-frame "$first_frame" "$@"
fi
critical sed -i -e '/nan/d' -e '/inf/d' "$output"
mark_done "${name}_density_analysis${suffix}"

