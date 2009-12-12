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
This script make all the post update with backup for single pairs

Usage: ${0##*/} step_nr

USES:  csg_get_interaction_property log mv die cp do_external run_or_exit

NEEDS: name inverse.post_update
EOF
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

name=$(csg_get_interaction_property name)
tasklist=$(csg_get_interaction_property --allow-empty inverse.post_update) 
i=1
for task in $tasklist; do
  log "Doing $task for ${name}"
  run_or_exit mv ${name}.dpot.new ${name}.dpot.cur
  run_or_exit cp ${name}.dpot.cur ${name}.dpot.${i}
  do_external postupd "$task" "$1"
  ((i++))
done
