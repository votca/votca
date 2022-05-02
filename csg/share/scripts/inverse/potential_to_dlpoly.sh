#!/bin/bash
#
# Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
trunc=${1##*/}
trunc="${trunc%%.*}"
[[ -f $input ]] || die "${0##*/}: Could not find input file '$input'"
output="$2"
echo "Convert $input to $output"

table_zero="0"
step="$(csg_get_interaction_property step)"
bondtype="$(csg_get_interaction_property bondtype)"
if [[ $bondtype = "non-bonded" ]]; then 
  OUT="TABLE"
  table_end="$(csg_get_property cg.inverse.dlpoly.table_end)"
  table_grid="$(csg_get_property cg.inverse.dlpoly.table_grid)"
  bin_size="$(csg_calc "$table_end" "/" $table_grid)"
  # make sure the TAB* file is removed externally
  if [[ ! -f "${OUT}" ]]; then
    echo "Table for dlpoly from VOTCA with love" > "${OUT}" #max 80 chars
    # see dlpoly manual ngrid = int(cut/delta) + 4
    table_grid="$(($table_grid+4))"
    # nm -> Angs
    bin_size1="$(csg_calc "$bin_size" "*" 10)"
    table_end1="$(csg_calc "$table_end" "*" 10)"
    echo "$bin_size1 $table_end1 $table_grid" >> "${OUT}"
  fi
elif [[ $bondtype = "bond" ]]; then 
  OUT="TABBND"
  table_end="$(csg_get_property cg.inverse.dlpoly.bonds.table_end)"
  table_grid="$(csg_get_property cg.inverse.dlpoly.bonds.table_grid)"
  bin_size="$(csg_calc "$table_end" "/" $table_grid)"
  # make sure the TAB* file is removed externally
  if [[ ! -f "${OUT}" ]]; then
    echo "# Table for dlpoly from VOTCA with love" > "${OUT}" #max 80 chars
    # nm -> Angs
    table_end1="$(csg_calc "$table_end" "*" 10)"
    echo "# $table_end1 $table_grid" >> "${OUT}"
  fi
elif [[ $bondtype = "angle" ]]; then 
  OUT="TABANG"
  table_end="3.14159265359"
  table_grid="$(csg_get_property cg.inverse.dlpoly.angles.table_grid)"
  bin_size="$(csg_calc "$table_end" "/" $table_grid)"
  # make sure the TAB* file is removed externally
  if [[ ! -f "${OUT}" ]]; then
    echo "# Table for dlpoly from VOTCA with love" > "${OUT}" #max 80 chars
    echo "# $table_grid" >> "${OUT}"
  fi
elif [[ $bondtype = "dihedral" ]]; then
  OUT="TABDIH"
  table_zero="-3.14159265359"
  table_end="3.14159265359"
  table_grid="$(csg_get_property cg.inverse.dlpoly.dihedrals.table_grid)"
  bin_size="$(csg_calc "$table_end" "-" $table_zero)"
  bin_size="$(csg_calc "$bin_size" "/" $table_grid)"
  # make sure the TAB* file is removed externally
  if [[ ! -f "${OUT}" ]]; then
    echo "# Table for dlpoly from VOTCA with love" > "${OUT}" #max 80 chars
    echo "# $table_grid" >> "${OUT}"
  fi
else
  die "${0##*/}: conversion of ${bondtype} interaction to generic tables is not implemented yet!"
fi
# Yes, the dlpoly table starts at ${bin_size}
table_begin="$(csg_calc "$table_zero" "+" $bin_size)"

#keep the grid for now, so that extrapolate can calculate the right mean
comment="$(get_table_comment)"
smooth2="$(critical mktemp ${trunc}.pot.extended.XXXXX)"
critical csg_resample --in ${input} --out "${smooth2}" --grid "${table_zero}:${step}:${table_end}" --comment "$comment"
extrapolate="$(critical mktemp ${trunc}.pot.extrapolated.XXXXX)"
do_external potential extrapolate --type "$bondtype" "${smooth2}" "${extrapolate}"

smooth="$(critical mktemp ${trunc}.pot.smooth.XXXXX)"
deriv="$(critical mktemp ${trunc}.pot.deriv.XXXXX)"
critical csg_resample --in ${extrapolate} --out "${smooth}" --der "${deriv}" --grid "${table_begin}:${bin_size}:${table_end}" --comment "$comment"

#shift does not change derivative
tshift="$(critical mktemp ${trunc}.pot.shift.XXXXX)"
do_external potential shift --type "$bondtype" "${smooth}" "${tshift}"

do_external convert_potential tab --header dlpoly --type "${bondtype}" "${tshift}" "${deriv}" "${output}"

if [[ -f $OUT ]]; then
  echo "Appending $output to $OUT" 
  if [[ $bondtype = "non-bonded" ]]; then
    #votca non-bonded types might not correspond to dl_poly's internal types, only use a failback
    header="$(csg_get_interaction_property --allow-empty dlpoly.header)"
    [[ -z ${header} ]] && header="$(csg_get_interaction_property type1) $(csg_get_interaction_property type2)"
    echo "${header}" >> "$OUT"
  else
    header="$(csg_get_interaction_property dlpoly.header)"
    # an empty line must precede each data block (for another bond type), then the bond type (two atom types) follow
    echo "" >> "$OUT"
    echo "# ${header}" >> "$OUT"
  fi
  cat "${output}" >> "$OUT"
fi
