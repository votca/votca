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
postadd overwrite script, overwrites potential of all other interactions with this one

Usage: ${0##*/} infile outfile
EOF
   exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"
critical cp "$1" "$2"

name=$(csg_get_interaction_property name)
tasklist=$(csg_get_interaction_property inverse.post_add)


#actually not possible, still checking
[[ -z $tasklist ]] && die "${0##*/}: Strange - how do you get here ?"

#has to be called at the end of postadd tasks and only at the end
[[ $tasklist =~ overwrite[[:space:]]*$ ]] || die "${0##*/}: Overwrite has to be called the last action in inverse.post_add"

step_nr=$(get_current_step_nr)
scheme=( $(csg_get_interaction_property inverse.post_add_options.overwrite.do 1) )
scheme_nr=$(( ($step_nr - 1 ) % ${#scheme[@]} ))

#this interaction will overwrite the others
if [ "${scheme[$scheme_nr]}" = 1 ]; then
  names="$(csg_get_property cg.non-bonded.name)"
  for i in $names; do
    [[ $i = $name ]] && continue
    echo "Overwriting $i.pot.new with $1"
    critical cp "$1" "$i.pot.new"
    mark_done "post_add-$i"
  done
fi
