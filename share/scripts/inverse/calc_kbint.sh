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

clean=no

show_help () {
  cat <<EOF
${0##*/}, version %version%
This script calculates the Kirkwood-Buff integral out of the rdf 

Usage: ${0##*/} [options] infile outfile
Allowed options:
    --help                    show this help
    --clean                   remove all intermediate temp files
EOF
}

### begin parsing options
shopt -s extglob
while [[ ${1} = --* ]]; do
 case $1 in
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

[[ -f $1 ]] || die "${0##*/}: $1 not found"

tmpfile=$(critical mktemp "$1".preint.XXX)
four_pi="12.5663706"

do_external table linearop "$1" "$tmpfile" "${four_pi}" "-${four_pi}"
do_external table integrate --from left --sphere "$tmpfile" "$2"

if [[ $clean = "yes" ]]; then
  rm -f "${tmpfile}"
fi
