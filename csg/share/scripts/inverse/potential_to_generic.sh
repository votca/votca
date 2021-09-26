#!/bin/bash
#
# Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
This script is a high class wrapper to convert a potential to the generic
3 column tab format used by espresso and espressopp

Usage: ${0##*/}
EOF
  exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: missing argument"
input="$1"
trunc=${1##*/}
trunc="${trunc%%.*}"
[[ -f $input ]] || die "${0##*/}: Could not find input file '$input'"
output="$2"
echo "Convert $input to $output"

bondtype="$(csg_get_interaction_property bondtype)"
[[ $bondtype != non-bonded ]] && die "${0##*/}: conversion of bonded interaction to generic tables is not implemented yet!"

sim_prog="$(csg_get_property cg.inverse.program)"
r_cut=$(csg_get_interaction_property max)
step=$(csg_get_interaction_property step)
r_min=$(csg_get_interaction_property min)
bin_size="$(csg_get_interaction_property --allow-empty inverse.$sim_prog.table_bins)"
[[ -z ${bin_size} ]] && bin_size="${step}"

comment="$(get_table_comment)"
table_begin="$(csg_get_interaction_property --allow-empty inverse.$sim_prog.table_begin)"
table_end="$(csg_get_interaction_property --allow-empty inverse.$sim_prog.table_end)"
[[ -z ${table_end} ]] && table_end="${r_cut}"
if [[ -z ${table_begin} ]]; then 
  table_begin="$r_min"
else
  #keep the grid for now, so that extrapolate can calculate the right mean
  smooth2="$(critical mktemp ${trunc}.pot.extended.XXXXX)"
  critical csg_resample --in ${input} --out "${smooth2}" --grid "${table_begin}:${step}:${table_end}" --comment "$comment"
  extrapolate="$(critical mktemp ${trunc}.pot.extrapolated.XXXXX)"
  lfct="$(csg_get_interaction_property --allow-empty inverse.$sim_prog.table_left_extrapolation)"
  rfct="$(csg_get_interaction_property --allow-empty inverse.$sim_prog.table_right_extrapolation)"
  do_external potential extrapolate ${lfct:+--lfct ${lfct}} ${rfct:+--rfct ${rfct}} --type "$bondtype" "${smooth2}" "${extrapolate}"
  input="${extrapolate}"
fi

#enable shift whenever bonded interaction are supported here
#do_external potential shift --type "$bondtype" "${extrapolate}" "${tshift}"

smooth="$(critical mktemp ${trunc}.pot.smooth.XXXXX)"
deriv="$(critical mktemp ${trunc}.pot.deriv.XXXXX)"
critical csg_resample --in ${input} --out "${smooth}" --der "${deriv}" --grid "${table_begin}:${bin_size}:${table_end}" --comment "$comment"
do_external convert_potential tab --header "${sim_prog}" --type "${bondtype}" "${smooth}" "${deriv}" "${output}"
