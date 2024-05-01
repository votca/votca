#! /usr/bin/env bash
#
# Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

clean=no
y1=0
y2=0

show_help () {
  cat <<EOF
${0##*/}, version %version%
This script creates a zero table with grid min:step:max using linear interpolation

Usage: ${0##*/} [options] min:step:max outfile
Allowed options:
    --y1  X.X                 using X.X instead of 0 for the 1st y-value
                              this creates a linear instead of a constant table
    --y2  X.X                 using X.X instead of 0 for the 2nd y-value
                              this creates a linear instead of a constant table
    --help                    show this help
    --clean                   remove all intermediate temp files
EOF
}

### begin parsing options
shopt -s extglob
while [[ ${1} = --* ]]; do
 case $1 in
   --y1)
    is_num "$2" || die "${0##*/}: argument of --y1 should be a number, but found '$2'"
    y1="$2"
    shift 2;;
   --y2)
    is_num "$2" || die "${0##*/}: argument of --y2 should be a number, but found '$2'"
    y2="$2"
    shift 2;;
   --clean)
    clean="yes"
    shift ;;
  --help)
    show_help
    exit 0;;
  *)
   die "Unknown option '$1'";;
 esac
done
[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"

if [[ ${1//[^:]} = "::" ]]; then
  min=${1%%:*}
  max=${1##*:}
else
  die "${0##*/}: Agrument 1 should have the form XX:XX:XX"
fi

tmpfile="$(critical mktemp table.XXX)"
echo "$min ${y1}" > $tmpfile
echo "$max ${y2}" >> $tmpfile

comment="$(get_table_comment)"
critical csg_resample --type linear --in ${tmpfile} --out "${2}" --grid "${1}" --comment "${comment}"

if [[ $clean = "yes" ]]; then
  rm -f "${tmpfile}"
fi
