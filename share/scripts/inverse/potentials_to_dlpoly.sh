#! /bin/bash
#
# Copyright 2009-2014 The VOTCA Development Team (http://www.votca.org)
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
This script converts all potentials to the format needed by dlpoly

Usage: ${0##*/}
EOF
   exit 0
fi

for i in TABLE TABBND TABANG TABDIH; do
  if [[ -f "${i}" ]]; then
    echo "We will now overwrite ${i}"
    rm -v "${i}"
  fi
done

for_all "non-bonded" touch "TABLE"
for_all "bond"       touch "TABBND"
for_all "angle"      touch "TABANG"
for_all "dihedral"   touch "TABDIH"

#if we have at least one  interaction for that kind
[[ -f "TABLE" ]] && echo "Table for dlpoly from VOTCA with love" > "TABLE" #max 80 chars
for i in TABBND TABANG TABDIH; do
  [[ -f "${i}" ]] && echo "# Table for dlpoly from VOTCA with love" > "${i}" #max 80 chars
done

if [[ -f "TABLE" ]]; then
  bin_size="$(csg_get_property cg.inverse.dlpoly.table_bins)"
  table_end="$(csg_get_property cg.inverse.dlpoly.table_end)"

  # see dlpoly manual ngrid = int(cut/delta) + 4
  ngrid="$(csg_calc $table_end / $bin_size)"
  ngrid="$(to_int $ngrid)"
  ngrid="$(($ngrid+4))"

  # nm -> Angs
  bin_size="$(csg_calc "$bin_size" "*" 10)"
  table_end="$(csg_calc "$table_end" "*" 10)"
  echo "$bin_size $table_end $ngrid" >> "TABLE"
  for_all "non-bonded" do_external convert_potential dlpoly '$(csg_get_interaction_property name).pot.cur' '$(csg_get_interaction_property name).pot.dlpoly'
fi

if [[ -f "TABBND" ]]; then
  table_end="$(csg_get_property cg.inverse.dlpoly.bonds.table_end)"
  ngrid="$(csg_get_property cg.inverse.dlpoly.bonds.table_grid)"

  # nm -> Angs
  table_end="$(csg_calc "$table_end" "*" 10)"
  echo "# $table_end $ngrid" >> "TABBND"
  for_all "bond" do_external convert_potential dlpoly '$(csg_get_interaction_property name).pot.cur' '$(csg_get_interaction_property name).pot.dlpoly'
fi

if [[ -f "TABANG" ]]; then
  ngrid="$(csg_get_property cg.inverse.dlpoly.angles.table_grid)"
  echo "# $ngrid" >> "TABANG"
  for_all "angle" do_external convert_potential dlpoly '$(csg_get_interaction_property name).pot.cur' '$(csg_get_interaction_property name).pot.dlpoly'
fi

if [[ -f "TABDIH" ]]; then
  ngrid="$(csg_get_property cg.inverse.dlpoly.dihedrals.table_grid)"
  echo "# $ngrid" >> "TABDIH"
  for_all "dihedral" do_external convert_potential dlpoly '$(csg_get_interaction_property name).pot.cur' '$(csg_get_interaction_property name).pot.dlpoly'
fi
