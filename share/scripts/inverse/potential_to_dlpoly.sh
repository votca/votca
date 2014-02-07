#!/bin/bash
#
# Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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
This script is a high class wrapper to convert a potential to the dlpoly format

Usage: ${0##*/}
EOF
  exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: missing argument"
input="$1"
trunc="${1%%.*}"
[[ -f $input ]] || die "${0##*/}: Could not find input file '$input'"
output="$2"
echo "Convert $input to $output"

bondtype="$(csg_get_interaction_property bondtype)"
[[ $bondtype != non-bonded ]] && die "${0##*/}: conversion of bonded interaction to generic tables is not implemented yet!"

step=$(csg_get_interaction_property step)
bin_size="$(csg_get_property cg.inverse.dlpoly.table_bins)"
table_end="$(csg_get_property cg.inverse.dlpoly.table_end)"

#keep the grid for now, so that extrapolate can calculate the right mean
comment="$(get_table_comment)"
smooth2="$(critical mktemp ${trunc}.pot.extended.XXXXX)"
critical csg_resample --in ${input} --out "${smooth2}" --grid "0:${step}:${table_end}" --comment "$comment"
extrapolate="$(critical mktemp ${trunc}.pot.extrapolated.XXXXX)"
do_external potential extrapolate --type "$bondtype" "${smooth2}" "${extrapolate}"

smooth="$(critical mktemp ${trunc}.pot.smooth.XXXXX)"
deriv="$(critical mktemp ${trunc}.pot.deriv.XXXXX)"
# Yes, the dlpoly table starts at ${bin_size}
critical csg_resample --in ${extrapolate} --out "${smooth}" --der "${deriv}" --grid "${bin_size}:${bin_size}:${table_end}" --comment "$comment"
do_external convert_potential tab --header dlpoly --type "${bondtype}" "${smooth}" "${deriv}" "${output}"

OUT="TABLE"
[[ $bondtype = bond ]] && OUT="TABBND"
[[ $bondtype = angle ]] && OUT="TABANG"
[[ $bondtype = dihedral ]] && OUT="TABDIH"

if [[ -f $OUT ]]; then
  echo "Appending $output to $OUT"
  if [[ $bondtype = non-bonded ]]; then
    #votca types might not correspond to dl_poly's internal types
    header="$(csg_get_interaction_property --allow-empty dlpoly.header)"
    [[ -z ${header} ]] && header="$(csg_get_interaction_property type1) $(csg_get_interaction_property type2)"
    echo "${header}" >> "$OUT"
  else
    echo "$(csg_get_interaction_property dlpoly.header)" >> "$OUT"
  cat "${output}" >> "$OUT"
fi
