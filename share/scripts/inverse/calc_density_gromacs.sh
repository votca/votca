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
max=$(csg_get_interaction_property tf.spline_end)
step=$(csg_get_interaction_property step)
bins=$(csg_calc $max / $step )

adress_type=$(get_simulation_setting adress_type)
echo "Adress type: $adress_type"

equi_time="$(csg_get_property cg.inverse.$sim_prog.equi_time 0)"
first_frame="$(csg_get_property cg.inverse.$sim_prog.first_frame 0)"
mol="$(csg_get_interaction_property tf.molname "*")"
if [ "$adress_type" = "sphere" ]; then
  adressc="$(get_simulation_setting adress_reference_coords "0 0 0")"
  ref="$(echo "$adressc" | awk '{if (NF<3) exit 1; printf "[%s,%s,%s]",$1,$2,$3;}')" || die "${0##*/}: we need three numbers in adress_reference_coords, but got '$adressc'"
  axis="r"
  opts="--ref $ref"
else
  axis="x"
fi

msg "Calculating density for $name (molname $mol) on axis $axis"
if is_done "density_analysis-$name"; then
  echo "density analysis is already done"
else
  critical csg_density --trj "$traj" --top "$topol" --out "$name.dist.new" --begin "$equi_time" --first-frame "$first_frame" --rmax "$max" --bins "$bins" --axis "$axis" --molname "$mol" $opts
  critical sed -i -e '/nan/d' -e '/inf/d' "$name.dist.new"
  mark_done "density_analysis-$name"
fi
