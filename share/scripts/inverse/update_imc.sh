#! /bin/bash
#
# Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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
This script implements the function update
for the Inverse Monte Carlo Method

Usage: ${0##*/}
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"
do_external imc_stat $sim_prog

default_reg=$(csg_get_property cg.inverse.imc.default_reg)
is_num "${default_reg}" || die "${0##*/}: value of cg.inverse.imc.default_reg should be a number"

imc_groups=$(csg_get_interaction_property --all inverse.imc.group)
imc_groups=$(remove_duplicate $imc_groups)
[[ -z ${imc_groups} ]] && die "${0##*/}: No imc groups defined"
for group in $imc_groups; do
  reg="$(csg_get_property cg.inverse.imc.${group}.reg ${default_reg})" #filter me away
  is_num "${reg}" || die "${0##*/}: value of cg.inverse.imc.${group}.reg should be a number"
  msg "solving linear equations for imc group '$group' (regularization ${reg})"
  critical csg_imc_solve --imcfile "${group}.imc" --gmcfile "${group}.gmc" --idxfile "${group}.idx" --regularization "${reg}" --outputfile "${group}.sol"
done

for_all "non-bonded bonded" do_external imc purify
