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
This script makes all the post update with backup for single pairs

Usage: ${0##*/}
EOF
   exit 0
fi

name=$(csg_get_interaction_property name)

#could be done by a overwrite somewhere else
is_done "post_add-$name" && exit 0

tasklist=$(csg_get_interaction_property --allow-empty inverse.post_add)
[[ -n $tasklist ]] && msg "Postadd tasks for $name: $tasklist"
i=1
#after all we tag the potential
for task in tag $tasklist; do
  echo "Doing postadd task '$task' for '${name}'"
  
  #save the current one
  critical mv "${name}.pot.new" "${name}.pot.${i}"

  #perform postadd task
  do_external postadd "$task" "${name}.pot.${i}" "${name}.pot.new" 

  ((i++))
done
mark_done "post_add-$name"
