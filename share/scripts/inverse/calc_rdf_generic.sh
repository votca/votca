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
This script implemtents statistical analysis for the iterative Boltzmann inversion 
using generic csg tools (csg_stat)

Usage: ${0##*/}
EOF
   exit 0
fi

name="$(csg_get_interaction_property name)"
sim_prog="$(csg_get_property cg.inverse.program)"

if [[ $sim_prog = "gromacs" ]]; then
  topol=$(csg_get_property cg.inverse.gromacs.topol_out)
  topol=$(csg_get_property cg.inverse.gromacs.rdf.topol "$topol")
  [[ -f $topol ]] || die "${0##*/}: gromacs topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist" 
  ext=$(csg_get_property cg.inverse.gromacs.traj_type)
  traj="traj.${ext}"
  [[ -f $traj ]] || die "${0##*/}: gromacs traj file '$traj' not found"
else
  die "${0##*/}: Simulation program '$sim_prog' not supported yet"
fi

if [[ -n $(csg_get_property --allow-empty cg.bonded.name) ]]; then
  mapping="$(csg_get_property cg.inverse.map)"
  mapping="$(csg_get_property cg.inverse.gromacs.rdf.map "$mapping")"
  [[ -f "$(get_main_dir)/$mapping" ]] || die "${0##*/}: Mapping file '$mapping' for bonded interaction not found in maindir"
  mapping="--cg $(get_main_dir)/$mapping"
else
  mapping=""
fi

equi_time="$(csg_get_property cg.inverse.gromacs.equi_time)"
first_frame="$(csg_get_property cg.inverse.gromacs.first_frame)"

with_errors=$(csg_get_property cg.inverse.gromacs.rdf.with_errors)
if [[ ${with_errors} = "yes" ]]; then
  suffix="_with_errors"
  block_length=$(csg_get_property cg.inverse.gromacs.rdf.block_length)
  error_opts="--block-length ${block_length} --ext dist.block"
else
  suffix=""
fi

tasks=$(get_number_tasks)
#rdf calculation is maybe done already in a different interaction
if is_done "rdf_calculation${suffix}"; then
  echo "rdf calculation is already done"
else
  msg "Calculating rdfs with csg_stat using $tasks tasks"
  critical csg_stat --nt $tasks --options "$CSGXMLFILE" --top "$topol" --trj "$traj" --begin $equi_time --first-frame $first_frame ${error_opts} ${mapping}
  mark_done "rdf_calculation${suffix}"
fi

if [[ ${with_errors} = "yes" ]]; then
  if ! is_done "${name}_rdf_average"; then
    msg "Calculating average rdfs and its errors for interaction $name"
    do_external table average --output ${name}.dist.new ${name}_*.dist.block
    mark_done "${name}_rdf_average"
  fi
fi
