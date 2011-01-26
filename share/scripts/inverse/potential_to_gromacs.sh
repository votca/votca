#!/bin/bash
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
This script is a wrapper to convert a potential to gromacs

Usage: ${0##*/}

USES: do_external csg_get_interaction_property csg_get_property critical csg_resample check_deps get_from_mdp awk

EOF
  exit 0
fi

check_deps "$0"

name=$(csg_get_interaction_property name)
input="${name}.pot.cur"
#gromacs want '_' !
output="$(csg_get_interaction_property inverse.gromacs.table)"
echo "Convert $input to $output"

r_cut=$(csg_get_interaction_property max)
gromacs_bins="$(csg_get_property cg.inverse.gromacs.table_bins)"

comment="$(get_table_comment)"

mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
[ -f "$mdp" ] || die "${0##*/}: gromacs mdp file '$mdp' not found"

rlist=$(get_from_mdp rlist "$mdp" 1)
tabext=$(get_from_mdp table-extension "$mdp" 1)

tablend="$(csg_get_property cg.inverse.gromacs.table_end)"
res="$(awk -v rlist="$rlist" -v tabext="$tabext" -v tablend="$tablend" 'BEGIN{ print (tablend<rlist+tabext)?1:0 }')" || die "${0##*/}: awk failed"
[ "$res" = "0" ] || die "${0##*/}: Error table is short then what gromacs needs, increase cg.inverse.gromacs.table_end in setting file.\nrlist ($rlist) + tabext ($tabext) > cg.inverse.gromacs.table_end ($tablend)"

critical csg_resample --in ${input} --out smooth_${input} --grid 0:${gromacs_bins}:${r_cut} --comment "$comment"
do_external convert_potential xvg smooth_${input} ${output}
