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
This script implemtents smoothing of the potential update (.dpot)

Usage: ${0##*/} infile outfile

USES: die csg_get_interaction_property mktemp sed awk csg_resample check_deps

NEEDS: name min max step inverse.post_update_options.splinesmooth.step
EOF
   exit 0
fi

check_deps "$0"

[ -z "$2" ] && die "${0##*/}: Missing arguments"

[ -f "$2" ] && die "${0##*/}: $2 is already there"

name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)
step=$(csg_get_interaction_property step)

tmpfile=$(mktemp ${name}.XXX) || die "mktemp failed"

sed -ne '/i[[:space:]]*$/p' "$1" > $tmpfile
spmin=$(sed -ne '1p' $tmpfile | awk '{print $1}')
spmax=$(sed -ne '$p' $tmpfile | awk '{print $1}')
spstep=$(csg_get_interaction_property inverse.post_update_options.splinesmooth.step)

comment="$(get_table_comment)"
critical csg_resample --in $tmpfile --out "$2" --grid $min:$step:$max --spfit $spmin:$spstep:$spmax --comment "$comment"

