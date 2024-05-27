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
This script initializes potentials for imc

Usage: ${0##*/}
EOF
   exit 0
fi

names=( $(csg_get_interaction_property --all name) )
if [[ ${#names[@]} -gt 1 ]]; then
  msg --color blue "####################################################"
  msg --color blue "# WARNING multicomponent imc is still experimental #"
  msg --color blue "####################################################"
fi

check_bonded_update() {
  local do_potential=$(csg_get_interaction_property inverse.do_potential)
  local imc_group=$(csg_get_interaction_property inverse.imc.group)
  if [[ $do_potential != 0 ]] || [[ $imc_group != "none" ]]; then
    die_msg="using IMC for bonded potentials is not implemented yet.\n"\
"Make sure to set update_potential to 0 and imc.group to 'none' for\n"\
"each bonded interaction. You can still use post_update ibi."
    die "$die_msg"
  fi
}
export -f check_bonded_update
for_all "bonded" check_bonded_update

do_external prepare generic
