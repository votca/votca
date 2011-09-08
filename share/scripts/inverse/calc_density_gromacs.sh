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
This script calcs the density for gromacs for the AdResS therm force

Usage: ${0##*/}
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"

if [ "$sim_prog" = "gromacs" ]; then
  topol=$(csg_get_property cg.inverse.gromacs.topol "topol.tpr")
  [ -f "$topol" ] || die "${0##*/}: gromacs topol file '$topol' not found"

  ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
  traj="traj.${ext}"
  [ -f "$traj" ] || die "${0##*/}: gromacs traj file '$traj' not found"
else
  die "${0##*/}: Simulation program '$sim_prog' not supported yet"
fi

name=$(csg_get_interaction_property name)

adress_type=$(get_simulation_setting adress_type)

equi_time="$(csg_get_property cg.inverse.$sim_prog.equi_time 0)"
first_frame="$(csg_get_property cg.inverse.$sim_prog.first_frame 0)"
mol="$(csg_get_interaction_property tf.molname "*")"
if [[ $adress_type = "sphere" ]]; then
  echo "Adress type: $adress_type"
  max=$(csg_get_interaction_property tf.spline_end)
  step=$(csg_get_interaction_property step)
  bins=$(csg_calc $max / $step )
  adressc="$(get_simulation_setting adress_reference_coords "0 0 0")"
  ref="$(echo "$adressc" | awk '{if (NF<3) exit 1; printf "[%s,%s,%s]",$1,$2,$3;}')" || die "${0##*/}: we need three numbers in adress_reference_coords, but got '$adressc'"
  axis="r"
  opts="--rmax $max --ref $ref"
#elif [[ $adress_type = "Xsplit" ]]
  else
  echo "Adress type: $adress_type"
  axis="x"
  opts=""
  bins="$(csg_get_interaction_property tf.density.bins)"
fi

with_errors=$(csg_get_property cg.inverse.gromacs.density.with_errors "no")
if [[ ${with_errors} = "yes" ]]; then
  suffix="_with_errors"
  output="$name.dist.block"
else
  suffix=""
  output="$name.dist.new"
fi

if is_done "${name}_density_analysis${suffix}"; then
  echo "density analysis is already done"
  exit 0
fi

with_errors=$(csg_get_property cg.inverse.gromacs.density.with_errors "no")
if [[ ${with_errors} = "yes" ]]; then
  msg "Calculating density for $name (molname $mol) on axis $axis with errors"
  block_length=$(csg_get_property cg.inverse.gromacs.density.block_length)
  critical csg_density --trj "$traj" --top "$topol" --out "$name.dist.block" --begin "$equi_time" --first-frame "$first_frame" --bins "$bins" --axis "$axis" --molname "$mol" --block-length $block_length $opts
  #mind the --clean option to avoid ${name}.dist.block_* to fail on the second run
  do_external table average --clean --output ${name}.dist.new ${name}.dist.block_*
else
  msg "Calculating density for $name (molname $mol) on axis $axis"
  critical csg_density --trj "$traj" --top "$topol" --out "$name.dist.new" --begin "$equi_time" --first-frame "$first_frame" --bins "$bins" --axis "$axis" --molname "$mol" $opts
  critical sed -i -e '/nan/d' -e '/inf/d' "$name.dist.new"
fi
mark_done "${name}_density_analysis${suffix}"

