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
Add table_comment to the head of a file

Usage: ${0##*/} input output

USES:  die msg do_external check_deps get_table_comment sed cat

NEEDS:
EOF
   exit 0
fi

check_deps "$0"

[[ -n "$2" ]] || die "${0##*/}: Missing arguments"

input="$1"
[ -f "$input" ] || die "${0##*/}: Input file $input not found"
output="$2"
[ -f "$output" ] && die "${0##*/}: $output is already there"

comment="$(get_table_comment)"

echo "Taging file $input to $output"
echo -e "$comment" | sed 's/^/#/' > "$output" || die "${0##*/}: sed failed"
cat "$input" >> "$output" || die "${0##*/}: sed failed"

