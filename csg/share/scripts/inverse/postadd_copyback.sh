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
postadd copyback script, copies files back to the maindir

Usage: ${0##*/}
EOF
   exit 0
fi

filelist=$(csg_get_interaction_property inverse.post_add_options.copyback.filelist)
echo "${0##*/}: copy $filelist to $(get_main_dir)"
critical cp $filelist $(get_main_dir)
