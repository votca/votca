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
This script adds up the tables

Usage: ${0##*/}
EOF
   exit 0
fi

msg "Adding up potential"
for_all "bonded non-bonded" \
  do_external table add '$(csg_get_interaction_property name).pot.cur $(csg_get_interaction_property name).dpot.new $(csg_get_interaction_property name).pot.new'
