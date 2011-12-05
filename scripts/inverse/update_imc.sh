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
This script implements the function update
for the Inverse Monte Carlo Method

Usage: ${0##*/}
EOF
   exit 0
fi

solver=$(csg_get_property cg.inverse.imc.solver)
sim_prog="$(csg_get_property cg.inverse.program)"
do_external imc_stat $sim_prog

imc_groups=$(csg_get_property cg.non-bonded.inverse.imc.group)
imc_groups=$(echo "$imc_groups" | sort -u)
for group in $imc_groups; do
  # currently this is a hack! need to create combined array
  msg "solving linear equations for $group"
  critical csg_imcrepack --in ${group} --out ${group}.packed
  do_external imcsolver $solver ${group}.packed ${group}.packed.sol
  critical csg_imcrepack --in ${group}.packed --unpack ${group}.packed.sol
done

for_all non-bonded do_external imc purify
