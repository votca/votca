#! /usr/bin/env bash
#
# Copyright 2009-2014 The VOTCA Development Team (http://www.votca.org)
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
This script cleans up after a simulation step

Usage: ${0##*/}
EOF
    exit 0
fi

cleanlist="$(csg_get_property --allow-empty cg.inverse.cleanlist)"
if [[ -n ${cleanlist} ]]; then
  msg "Clean up files: $cleanlist"
  #no quote to allow globbing
  rm -f ${cleanlist}
fi
