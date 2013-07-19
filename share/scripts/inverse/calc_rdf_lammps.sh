#! /bin/bash
#
# Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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
This script converts lammps' rdf output to VOTCA's internal format

Usage: ${0##*/}
EOF
    exit 0
fi

rdf="$(csg_get_property cg.inverse.lammps.rdf)"
[[ -f $rdf ]] || die "${0##*/}: rdf file '$rdf' could not be found."

type1=$(csg_get_interaction_property type1)
type2=$(csg_get_interaction_property type2)
name=$(csg_get_interaction_property name)
binsize=$(csg_get_interaction_property step)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)

echo "Analyzing rdf for ${type1}-${type2}"
if is_done "rdf-$name"; then
    echo "rdf analsysis for ${type1}-${type2} is already done"
else
    tmp=$(critical mktemp ${name}.dist.raw.XXXX)
    critical awk '/^#/{next;}(NF==4){print $2,$3}' "$rdf" > "${tmp}"

    comment="$(get_table_comment)"
    critical csg_resample --in "${tmp}" --out "${name}".dist.new --grid "${min}:${binsize}:${max}" --comment "$comment"
    mark_done "rdf-$name"
fi
