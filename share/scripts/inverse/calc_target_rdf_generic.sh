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
This script calculated reference rdf using generic csg_stat

Usage: ${0##*/}
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
[[ $sim_prog = "gromacs" ]] || die "${0##*/}: Simulation program '$sim_prog' not supported yet"

topol=$(csg_get_property cg.inverse.gromacs.ref.topol)
[[ $topol = /* ]] || topol=$(get_main_dir)/$topol 
[[ -f $topol ]] ||  die "${0##*/}: Reference toplogy '$topol' not found"
traj="$(csg_get_property cg.inverse.gromacs.ref.traj)"
[[ $traj = /* ]] || traj=$(get_main_dir)/$traj 
[[ -f $traj ]] ||  die "${0##*/}: Reference trajectory '$traj' not found"
mapping="$(csg_get_property cg.inverse.gromacs.ref.mapping)"
# no mapping check as this could be 1.xml;2.xml
equi_time="$(csg_get_property cg.inverse.gromacs.ref.equi_time)"
first_frame="$(csg_get_property cg.inverse.gromacs.ref.first_frame)"

opts="$(csg_get_property --allow-empty cg.inverse.gromacs.ref.rdf.opts)"

tasks=$(get_number_tasks)
#refrdf calculation is maybe done already in a different interaction
if is_done "refrdf_calculation"; then
  echo "reference rdf calculation is already done"
else
  msg "Calculating reference rdfs with csg_stat using $tasks tasks"
  critical csg_stat --nt $tasks --options "$CSGXMLFILE" --top "$topol" --trj "$traj" --begin $equi_time --first-frame $first_frame --cg "${mapping}" --ext "dist.tgt" ${opts}
  mark_done "refrdf_calculation"
fi
