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
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
This script adds LJ 12-6 component to the CG potential

Usage: ${0##*/} infile outfile
EOF
   exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"

[[ -f $2 ]] && die "${0##*/}: $2 is already there"

step_nr="$(get_current_step_nr)"

name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)
step=$(csg_get_interaction_property step)

echo "Add LJ12-6 component for interaction ${name}"
do_external compute_lj 12_6 ${name}.lj

comment="$(get_table_comment ${name}.lj)"
tmpfile=$(critical mktemp ${name}.lj.XXX)
critical csg_resample --in ${name}.lj --out ${tmpfile} --grid $min:$step:$max --comment "$comment"
do_external table add "$1" ${tmpfile} "$2"
