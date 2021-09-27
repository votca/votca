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
This script implemtents statistical analysis for the iterative Boltzmann inversion
using generic csg tools (csg_stat)

Usage: ${0##*/}
EOF
   exit 0
fi

name="$(csg_get_interaction_property name)"
sim_prog="$(csg_get_property cg.inverse.program)"

topol=$(csg_get_property --allow-empty cg.inverse.$sim_prog.rdf.topol)
[[ -z $topol ]] && topol=$(csg_get_property cg.inverse.$sim_prog.topol)
[[ -f $topol ]] || die "${0##*/}: topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist"

traj=$(csg_get_property --allow-empty cg.inverse.$sim_prog.rdf.traj)
[[ -z $traj ]] && traj=$(csg_get_property cg.inverse.$sim_prog.traj)
[[ -f $traj ]] || die "${0##*/}: traj file '$traj' not found"

maps=
#always try to find mapping files
if : ; then
  mapping="$(csg_get_property --allow-empty cg.inverse.$sim_prog.rdf.map)"
  [[ -z $mapping ]] && mapping="$(csg_get_property --allow-empty cg.inverse.map)"
  #fail if we have bonded interaction, but no mapping file
  [[ -n $(csg_get_property --allow-empty cg.bonded.name) && -z $mapping ]] && die "Mapping file for bonded interaction needed"
  for map in ${mapping}; do
    [[ -f "$(get_main_dir)/$map" ]] || die "${0##*/}: Mapping file '$map' for bonded interaction not found in maindir"
    maps+="$(get_main_dir)/$map;"
  done
fi

equi_time="$(csg_get_property cg.inverse.$sim_prog.equi_time)"
if [[ ${CSG_RUNTEST} ]] && csg_calc "$equi_time" ">" "0"; then
  msg --color blue --to-stderr "Automatically setting equi_time to 0, because CSG_RUNTEST was set"
  equi_time=0
fi

first_frame="$(csg_get_property cg.inverse.$sim_prog.first_frame)"
if [[ ${CSG_RUNTEST} ]] && csg_calc "$first_frame" ">" "0"; then
  msg --color blue --to-stderr "Automatically setting first_frame to 0, because CSG_RUNTEST was set"
  first_frame=0
fi

with_errors=$(csg_get_property cg.inverse.$sim_prog.rdf.with_errors)
if [[ ${with_errors} = "yes" ]]; then
  suffix="_with_errors"
  block_length=$(csg_get_property cg.inverse.$sim_prog.rdf.block_length)
  if [[ ${CSG_RUNTEST} ]] && csg_calc "$block_length" ">" "2"; then
    msg --color blue --to-stderr "Automatically setting block_length to 2, because CSG_RUNTEST was set"
    block_length=2
  fi
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
  critical csg_stat --nt $tasks --options "$CSGXMLFILE" --top "$topol" --trj "$traj" --begin $equi_time --first-frame $first_frame ${error_opts} ${maps:+--cg ${maps}}
  mark_done "rdf_calculation${suffix}"
fi

if [[ ${with_errors} = "yes" ]]; then
  if ! is_done "${name}_rdf_average"; then
    for i in ${name}_*.dist.block; do
      [[ -f $i ]] || die "${0##*/}: Could not find ${name}_*.dist.block after running csg_sat, that usually means the blocksize (cg.inverse.$sim_prog.rdf.block_length) is too big."
    done
    msg "Calculating average rdfs and its errors for interaction $name"
    do_external table average --output ${name}.dist.new ${name}_*.dist.block
    mark_done "${name}_rdf_average"
  fi
fi
