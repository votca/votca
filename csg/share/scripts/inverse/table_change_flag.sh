#! /usr/bin/env bash
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
This script changes the flags (col 3) of a table

Usage: ${0##*/} input outfile
EOF
   exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"

[[ -f $1 ]] || die "${0##*/}: Cannot find input file $1"

critical sed -e 's/o$/i/' $1 > $2
