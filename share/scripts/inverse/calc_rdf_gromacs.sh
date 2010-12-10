#! /bin/bash
#
# Copyright 2009 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicale law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script calcs the rdf for gromacs
for the Inverse Boltzmann Method

Usage: ${0##*/}

USES: get_from_mdp csg_get_interaction_property csg_get_property awk successful_or_die csg_resample is_done mark_done msg get_number_tasks multi_g_rdf check_deps

NEEDS: type1 type2 name step min max

OPTIONAL: cg.inverse.gromacs.equi_time cg.inverse.gromacs.first_frame cg.inverse.mpi.tasks cg.inverse.gromacs.mdp cg.inverse.gromacs.g_rdf.topol cg.inverse.gromacs.g_rdf.index cg.inverse.gromacs.g_rdf.opts cg.inverse.gromacs.traj_type cg.inverse.gromacs.g_rdf.bin
EOF
   exit 0
fi

check_deps "$0"

mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
[ -f "$mdp" ] || die "${0##*/}: gromacs mdp file '$mdp' not found"
dt=$(get_from_mdp dt "$mdp")
equi_time="$(csg_get_property cg.inverse.gromacs.equi_time 0)"
steps=$(get_from_mdp nsteps "$mdp")
first_frame="$(csg_get_property cg.inverse.gromacs.first_frame 0)"

index="$(csg_get_property cg.inverse.gromacs.g_rdf.index "index.ndx")"
[ -f "$index" ] || die "${0##*/}: grompp index file '$index' not found"
tpr="$(csg_get_property cg.inverse.gromacs.g_rdf.topol "topol.tpr")"
[ -f "$tpr" ] || die "${0##*/}: Gromacs tpr file '$tpr' not found"

ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
traj="traj.${ext}"
[ -f "$traj" ] || die "${0##*/}: gromacs traj file '$traj' not found"

g_rdf="$(csg_get_property cg.inverse.gromacs.g_rdf.bin "g_rdf")"
[ -n "$(type -p $g_rdf)" ] || die "${0##*/}: g_rdf binary '$g_rdf' not found"

opts="$(csg_get_property --allow-empty cg.inverse.gromacs.g_rdf.opts)"

type1=$(csg_get_interaction_property type1)
type2=$(csg_get_interaction_property type2)
name=$(csg_get_interaction_property name)
binsize=$(csg_get_interaction_property step)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)

begin="$(awk -v dt=$dt -v frames=$first_frame -v eqtime=$equi_time 'BEGIN{print (eqtime > dt*frames ? eqtime : dt*frames) }')"
end="$(awk -v dt="$dt" -v steps="$steps" 'BEGIN{print dt*steps}')"

tasks=$(get_number_tasks)
echo "Running g_rdf for ${type1}-${type2} using $tasks tasks"
if is_done "rdf-$name"; then
  msg "g_rdf for ${type1}-${type2} is already done"
else
  if [ $tasks -gt 1 ]; then
    echo -e "${type1}\n${type2}" | successful_or_die multi_g_rdf --cmd ${g_rdf} -${tasks} -b ${begin} -e ${end} -n "$index" -o ${name}.dist.new.xvg --soutput ${name}.dist.new.NP.xvg -- -bin ${binsize}  -s "$tpr" -f "${traj}" ${opts} 
  else
    echo -e "${type1}\n${type2}" | successful_or_die ${g_rdf} -b ${begin} -n "$index" -bin ${binsize} -o ${name}.dist.new.xvg -s "$tpr" -f "${traj}" ${opts}
  fi
  #gromacs always append xvg
  comment="$(get_table_comment)"
  successful_or_die csg_resample --in ${name}.dist.new.xvg --out ${name}.dist.new --grid ${min}:${binsize}:${max} --comment "$comment"
  mark_done "rdf-$name"
fi
