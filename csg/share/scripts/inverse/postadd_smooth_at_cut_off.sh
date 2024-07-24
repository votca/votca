#! /usr/bin/env bash
#
# Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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
This script implements smoothing of the potential (.pot) at the cut off

Usage: ${0##*/} infile outfile
EOF
   exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"

[ -f "$2" ] && die "${0##*/}: $2 is already there"

name=$(csg_get_interaction_property name)
tmpfile=$(critical mktemp "${name}.XXX")
cut_off=$(csg_get_interaction_property max)

critical cp "$1" "${tmpfile}"
echo "smoothing near cut-off for interaction ${name}"

critical do_external table smooth_at_cut_off "${tmpfile}" "$2" --cut-off="${cut_off}"
critical rm -f "${tmpfile}"
