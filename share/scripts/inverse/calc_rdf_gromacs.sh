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
# Unless required by applicale law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script calcs the rdf for gromacs using g_rdf

Usage: ${0##*/}

Used external packages: gromacs
EOF
   exit 0
fi

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

grp1=$(csg_get_interaction_property --allow-empty gromacs.grp1)
if [ -z "$grp1" ]; then
  grp1=$(csg_get_interaction_property type1)
  msg --color blue --to-stderr "WARNING (in ${0##*/}): could not find energy group of type1 (interaction property gromacs.grp1) using bead type1 ($grp1) as failback !"
fi
grp2=$(csg_get_interaction_property --allow-empty gromacs.grp2)
if [ -z "$grp2" ]; then
  grp2=$(csg_get_interaction_property type2)
  msg --color blue --to-stderr "WARNING (in ${0##*/}): could not find energy group of type2 (interaction property gromacs.grp2) using bead type2 ($grp2) as failback !"
fi
name=$(csg_get_interaction_property name)
binsize=$(csg_get_interaction_property step)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)

begin="$(awk -v dt=$dt -v frames=$first_frame -v eqtime=$equi_time 'BEGIN{print (eqtime > dt*frames ? eqtime : dt*frames) }')"
end="$(awk -v dt="$dt" -v steps="$steps" 'BEGIN{print dt*steps}')"

tasks=$(get_number_tasks)
if is_done "rdf-$name"; then
  echo "g_rdf for ${grp1}-${grp2} is already done"
else
  msg "Running g_rdf for ${grp1}-${grp2} using $tasks tasks"
  if [ $tasks -gt 1 ]; then
    echo -e "${grp1}\n${grp2}" | critical multi_g_rdf --cmd ${g_rdf} -${tasks} -b ${begin} -e ${end} -n "$index" -o ${name}.dist.new.xvg --soutput ${name}.dist.new.NP.xvg -- -bin ${binsize}  -s "$tpr" -f "${traj}" ${opts}
  else
    echo -e "${grp1}\n${grp2}" | critical ${g_rdf} -b ${begin} -n "$index" -bin ${binsize} -o ${name}.dist.new.xvg -s "$tpr" -f "${traj}" ${opts}
  fi
  #gromacs always append xvg
  comment="$(get_table_comment)"
  critical csg_resample --in ${name}.dist.new.xvg --out ${name}.dist.new --grid ${min}:${binsize}:${max} --comment "$comment"
  mark_done "rdf-$name"
fi
