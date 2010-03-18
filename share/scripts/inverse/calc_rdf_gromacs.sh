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
# Unless required by applicable law or agreed to in writing, software
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

USES: get_from_mdp csg_get_interaction_property csg_get_property awk log run_or_exit g_rdf csg_resample is_done mark_done msg use_mpi multi_g_rdf check_deps

NEEDS: type1 type2 name step min max

OPTIONAL: cg.inverse.gromacs.equi_time cg.inverse.gromacs.first_frame cg.inverse.mpi.tasks
EOF
   exit 0
fi

check_deps "$0"

dt=$(get_from_mdp dt "grompp.mdp")
equi_time="$(csg_get_property cg.inverse.gromacs.equi_time 0)"
steps=$(get_from_mdp nsteps "grompp.mdp")
first_frame="$(csg_get_property cg.inverse.gromacs.first_frame 0)"

type1=$(csg_get_interaction_property type1)
type2=$(csg_get_interaction_property type2)
name=$(csg_get_interaction_property name)
binsize=$(csg_get_interaction_property step)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)

begin="$(awk -v dt=$dt -v frames=$first_frame -v eqtime=$equi_time 'BEGIN{print (eqtime > dt*frames ? eqtime : dt*frames) }')"
end="$(awk -v dt="$dt" -v steps="$steps" 'BEGIN{print dt*steps}')"

log "Running g_rdf for ${type1}-${type2}"
if is_done "rdf-$name"; then
  msg "g_rdf for ${type1}-${type2} is already done"
else
  if use_mpi; then
    tasks=$(csg_get_property cg.inverse.mpi.tasks)
    echo -e "${type1}\n${type2}" | run_or_exit multi_g_rdf -${tasks} -b ${begin} -e ${end} -n index.ndx -o ${name}.dist.new.xvg --soutput ${name}.dist.new.NP.xvg -- -noxvgr -bin ${binsize}  -s topol.tpr
  else
    echo -e "${type1}\n${type2}" | run_or_exit g_rdf -b ${begin} -noxvgr -n index.ndx -bin ${binsize} -o ${name}.dist.new.xvg -s topol.tpr
  fi
#gromacs always append xvg
  comment="$(get_table_comment)"
  run_or_exit csg_resample --in ${name}.dist.new.xvg --out ${name}.dist.new --grid ${min}:${binsize}:${max} --comment "$comment"
  mark_done "rdf-$name"
fi
