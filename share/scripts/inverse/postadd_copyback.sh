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
postadd copyback script, copies files back to the maindir,
use \${name} in filename as replacement for the interaction name

Usage: ${0##*/} infile outfile
EOF
   exit 0
fi

[ -z "$2" ] && die "${0##*/}: Missing arguments"

do_external postadd dummy "$1" "$2"

filelist=$(csg_get_interaction_property --allow-empty inverse.post_add_options.copyback.filelist)
name=$(csg_get_interaction_property name)

if [ -n "$filelist" ]; then
  echo "${0##*/}: copy $filelist to $(get_main_dir)"
  #we do eval because i can contain $name
  critical cp_from_to --from . --where $(get_main_dir) $(eval echo $filelist)
fi

