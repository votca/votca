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

This script generates potential (.pot.new) for all interactions out the first pending line in the input state file and flags this line active in output state

Usage: ${0##*/} inputstate outputstate
EOF
  exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: missing argument"
input="$1"
ouput="$2"
[[ -f $input ]] || die "${0##*/}: Could not read $input"

line="$(critical sed -n '/pending$/{=;q}' "$input")"
[[ -z $line ]] && die "${0##*/}: not could find a pending line in $input"
is_int "$line" || die "${0##*/}: Strange - $line should be a number"
critical sed "${line}s/pending$/active/" "$input" > "$output"

parameters="$(sed -n "${line}p" "$input")"
for_all non-bonded do_external simplex parameters_to_potential "$parameters" 
