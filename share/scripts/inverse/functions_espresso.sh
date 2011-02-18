#!/bin/bash
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
Useful functions for espresso

Used external packages: espresso
EOF
  exit 0
fi

#check performed when this file is sourced
esp_dir="$(csg_get_property --allow-empty cg.inverse.espresso.scriptdir)"
if [ -z "${esp_dir}" ]; then
  [ -z "$ESPRESSO_SCRIPTS" ] && die "${0##*/}: cg.inverse.espresso.scriptdir of the xml setting file was empty and ESPRESSO_SCRIPTS not set in the environment.\nEspresso needs this variable to find its scripts."
  [ -d "${ESPRESSO_SCRIPTS}" ] || die "${0##*/}: ESPRESSO_SCRIPTS ($ESPRESSO_SCRIPTS) is not a directory"
else
  export ESPRESSO_SCRIPTS="${esp_dir}"
  [ -d "${ESPRESSO_SCRIPTS}" ] || die "${0##*/}: cg.inverse.espresso.scriptdir ($ESPRESSO_SCRIPTS) is not a directory"
fi
unset esp_dir

