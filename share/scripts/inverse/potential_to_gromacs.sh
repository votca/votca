#!/bin/bash
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
This script is a wrapper to convert a potential to gromacs

Usage: ${0##*/}
EOF
  exit 0
fi

name=$(csg_get_interaction_property name)
input="${name}.pot.cur"
#gromacs want '_' !
output="$(csg_get_interaction_property inverse.gromacs.table)"
echo "Convert $input to $output"

r_cut=$(csg_get_interaction_property max)
gromacs_bins="$(csg_get_property cg.inverse.gromacs.table_bins)"

comment="$(get_table_comment)"

tablend="$(csg_get_property --allow-empty cg.inverse.gromacs.table_end)"
mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
if [ -f "${mdp}" ]; then
  rlist=$(get_from_mdp rlist "$mdp" 0)
  tabext=$(get_from_mdp table-extension "$mdp" 0)
else
  rlist=0
  tabext=0
fi

if [ -z "${tablend}" ]; then
  ( [ "${rlist}" = "0" ] || [ "${tabext}" = "0" ] ) && die "${0##*/}: cg.inverse.gromacs.table_end not given and rlist and table-extension not found in mdp file ($mdp)"
  tablend=$(csg_calc "$rlist" + "$tabext")
elif [ "${rlist}" != "0" ] && [ "${tabext}" != "0" ]; then
  csg_calc $tablend "<" "$rlist+$tabext" && \
    die "${0##*/}: Error table is shorter then what gromacs needs, increase cg.inverse.gromacs.table_end in setting file.\nrlist ($rlist) + tabext ($tabext) > cg.inverse.gromacs.table_end ($tablend)"
fi

smooth="$(critical mktemp smooth_${name}.XXXXX)"
csg_resample --in ${input} --out "$smooth" --grid 0:${gromacs_bins}:${tablend} --comment "$comment"
extrapol="$(critical mktemp extrapol_${name}.XXXXX)"
do_external table extrapolate --function constant --avgpoints 1 --region leftright "${smooth}" "${extrapol}"
tshift="$(critical mktemp shift_${name}.XXXXX)"
#will shift the potential to the last defined value
do_external pot shift_nb "${extrapol}" "${tshift}"

potmax="$(csg_get_property --allow-empty cg.inverse.gromacs.pot_max)"
[ -n "${potmax}" ] && potmax="--max ${potmax}"
tabtype="$(csg_get_interaction_property bondtype)"
[ "$(csg_get_property cg.inverse.method)" = "tf" ] && tabtype="thermforce"
do_external convert_potential xvg ${potmax} --type "${tabtype}" "${tshift}" "${output}" 
