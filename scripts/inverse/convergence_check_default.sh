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
Calculated the sum of all convergence files and create a file 'stop' if the sum is bigger than a given limit

Usage: ${0##*/}
EOF
   exit 0
fi

limit="$(csg_get_property cg.inverse.convergence_check.limit)"

sum=0
names="$(csg_get_property cg.non-bonded.name)"
found=0
for name in $names; do
  out="${name}.conv"
  [[ -f $out ]] || continue
  ((found++))
  val="$(<$out)"
  is_num "$val" || die "${0##*/}: Content of $i was not a number"
  sum=$(csg_calc "$sum" + "$val")
done
[[ $found -eq 0 ]] && die "${0##*/}: No convergence file found!\nHave you added convergence to the postadd list of at least one interaction?"

echo "Convergence sum was $sum, limit is $limit"
csg_calc "$sum" ">" "$limit" || touch 'stop'
