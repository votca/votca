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

Usage: ${0##*/}

USES:  csg_get_interaction_property mv do_external critical

NEEDS: name inverse.post_add
EOF
   exit 0
fi

name=$(csg_get_interaction_property name)
tasklist=$(csg_get_interaction_property --allow-empty inverse.post_add)
i=1
#after all we shift the potential to be 0 at the cutoff and tag with labels
for task in $tasklist shift tag; do
  echo "Doing postadd task '$task' for '${name}'"
  
  #save the current one
  critical mv "${name}.pot.new" "${name}.pot.${i}"

  #perform postadd task
  do_external postadd "$task" "${name}.pot.${i}" "${name}.pot.new" 

  ((i++))
done
