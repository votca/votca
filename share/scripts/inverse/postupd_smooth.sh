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
This script implements smoothing of the potential update (.dpot)

Usage: ${0##*/} infile outfile
EOF
   exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"

[ -f "$2" ] && die "${0##*/}: $2 is already there"

name=$(csg_get_interaction_property name)
tmpfile=$(critical mktemp ${name}.XXX)
iterations=$(csg_get_interaction_property inverse.post_update_options.smooth.iterations 1)

critical cp "$1" $tmpfile
echo "doing $iterations smoothing iterations for interaction $name"

for((i=0;i<$iterations;i++)); do
  do_external table smooth $tmpfile "$2"
  critical cp "$2" $tmpfile
done

