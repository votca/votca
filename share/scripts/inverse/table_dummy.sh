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
if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script creates a dummy table with grid min:step:max

Usage: ${0##*/} min:step:max outfile
EOF
   exit 0
fi

[ -z "$2" ] && die "${0##*/}: Missing arguments"

if [ "${1//[^:]}" = "::" ]; then
  min=${1%%:*}
  max=${1##*:}
else
  die "${0##*/}: Agrument 1 should have the form XX:XX:XX"
fi

tmpfile="$(critical mktemp table.XXX)"
echo "$min 0" > $tmpfile
echo "$max 0" >> $tmpfile

comment="$(get_table_comment)"
critical csg_resample --type linear --in ${tmpfile} --out "${2}" --grid "${1}" --comment "${comment}"
