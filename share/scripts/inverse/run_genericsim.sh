#! /bin/bash
#
# Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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
This script runs a generic simulation program

Usage: ${0##*/}
EOF
    exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
script="$(csg_get_property --allow-empty cg.inverse.$sim_prog.script)"
[[ -n $script && ! -f $script ]] && die "${0##*/}: $sim_prog script '$script' not found (make sure it is in cg.inverse.filelist)"

cmd="$(csg_get_property cg.inverse.$sim_prog.command)"
#no check for cmd, because cmd could maybe exist only computenodes

opts=$(csg_get_property --allow-empty cg.inverse.$sim_prog.opts)
#expand ${script} in there
opts="$(eval echo $opts)"

if [[ -n $CSGENDING ]]; then
  echo "${0##*/}: $sim_prog does not support wallclock time yet (go here and implement it). Per step wallclock time check is still performed!"
fi

if [[ ${CSG_MDRUN_STEPS} ]]; then
  if [[ ${sim_prog} = "lammps" ]]; then
    critical sed -i "/^run/s/[0-9][0-9]*/${CSG_MDRUN_STEPS}/" "$script"
    msg --color blue --to-stderr "Replace run STEPS in '$script' to be ${CSG_MDRUN_STEPS}"
  elif [[ ${sim_prog} = "espresso" ]]; then
    critical sed -i -e "/^steps_per_int/s/[0-9][0-9]*/${CSG_MDRUN_STEPS}/" \
                    -e '/^\(int\|eq\)_steps/s/[0-9][0-9]*/1/' "$script"
    msg --color blue --to-stderr "Replace steps_per_int in '$script' to be ${CSG_MDRUN_STEPS} and int_steps to be 1"
  else
    msg --color blue --to-stderr "Overwriting of nsteps for ${sim_prog} not supported yet"
  fi
fi

method="$(csg_get_property cg.inverse.method)"
shopt -s extglob
is_part "$method" "ibi imc optimizer re" || die "${0##*/}: ${sim_prog} does not support method $method yet!"

critical $cmd ${opts}

simulation_finish
