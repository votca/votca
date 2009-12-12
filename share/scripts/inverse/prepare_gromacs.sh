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
This script implemtents the function initialize
for the Inverse Boltzmann Method

Usage: ${0##*/} last_sim_dir

USES: die cp run_or_exit grompp

NEEDS:
EOF
  exit 0
fi

check_deps "$0"

[[ -z "$1" ]] && die "Missing argument for ${0##*/}"

cp ${1}/confout.gro ./conf.gro || die "${0##*/} cp ${1}/confout.gro ./conf.gro failed" 

run_or_exit grompp -n index.ndx 
