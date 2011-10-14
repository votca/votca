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

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
This script implemtents the function initialize

Usage: ${0##*/}

Used external packages: gromacs
EOF
  exit 0
fi

from="$(csg_get_property cg.inverse.initial_configuration)"
conf="$(csg_get_property cg.inverse.gromacs.conf)"
if [[ $from = "laststep" ]]; then
  confout="$(csg_get_property cg.inverse.gromacs.conf_out)"
  #avoid overwriting $confout
  cp_from_last_step --rename "${confout}" "${conf}"
elif [[ $from = "maindir" ]]; then
  cp_from_main_dir "$conf"
else
  die "${0##*/}: initial_configuration '$from' not implemented"
fi

#convert potential in format for sim_prog
for_all "non-bonded bonded" do_external convert_potential gromacs

check_temp || die "${0##*/}: check of tempertures failed"

if [[ $(csg_get_property cg.inverse.method) != "tf" ]]; then
  for_all "non-bonded" check_cutoff || die "${0##*/}: check of cutoff for non-bonded interactions failed"
fi
