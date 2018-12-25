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

show_help () {
  cat <<EOF
${0##*/}, version %version%
This script is a high class wrapper to convert a potential to the lammps format

Usage: ${0##*/} [options] input output

Allowed options:
    --help       show this help
    --clean      remove all intermediate temp files
    --no-r2d     do not converts rad to degree (scale x axis with 180/3.1415)
                 for angle interactions
                 Note: VOTCA calcs in rad, but lammps uses degrees for angle
    --no-shift   do not shift the potential
EOF
}

clean=
do_shift="yes"
r2d="57.2957795"
step=

### begin parsing options
shopt -s extglob
while [[ ${1#-} != $1 ]]; do
 if [[ ${1#--} = $1 && -n ${1:2} ]]; then
    #short opt with arguments here: o
    if [[ ${1#-[o]} != ${1} ]]; then
       set -- "${1:0:2}" "${1:2}" "${@:2}"
    else
       set -- "${1:0:2}" "-${1:2}" "${@:2}"
    fi
 fi
 case $1 in
   --no-r2d)
    r2d=1
    shift ;;
   --clean)
    clean="yes"
    shift ;;
   --no-shift)
    do_shift="no"
    shift ;;
   -h | --help)
    show_help
    exit 0;;
  *)
   die "Unknown option '$1'";;
 esac
done

### end parsing options
[[ -z $1 || -z $2 ]] && die "${0##*/}: missing argument"
input="$1"
trunc=${1##*/}
trunc="${trunc%%.*}"
[[ -f $input ]] || die "${0##*/}: Could not find input file '$input'"
output="$2"
echo "Convert $input to $output"

sim_prog="$(csg_get_property cg.inverse.program)"

bondtype="$(csg_get_interaction_property bondtype)"

r_cut=$(csg_get_interaction_property max)
step=$(csg_get_interaction_property step)
r_min=$(csg_get_interaction_property min)
bin_size="$(csg_get_interaction_property --allow-empty inverse.$sim_prog.table_bins)"
[[ -z ${bin_size} ]] && bin_size="${step}"

if [[ $bondtype = "angle" ]]; then
  table_begin=0
  table_end=180
else
  table_begin="$(csg_get_interaction_property --allow-empty inverse.$sim_prog.table_begin)"
  [[ -z ${table_begin} ]] && table_begin="${r_min}"
  table_end="$(csg_get_interaction_property --allow-empty inverse.$sim_prog.table_end)"
  [[ -z ${table_end} ]] && table_end="${r_cut}"
fi

comment="$(get_table_comment $input)"

scale_factor=$(csg_get_interaction_property inverse.lammps.scale)

if [[ $bondtype = "angle" ]] && [[ $r2d != 1 ]]; then
# see lammps manual; only tabulated angle potentials need to be in degrees and scaled; dihedrals can be in radians or degrees; here we are putting dihedrals in radians
  scale="$(critical mktemp ${trunc}.pot.scale.XXXXX)"
  do_external table linearop --on-x "${input}" "${scale}" "$r2d" "0"
  step=$(csg_calc $r2d "*" $step)
elif [[ ${scale_factor} != 1 ]]; then
  step=$(csg_calc ${scale_factor} "*" $step)
  bin_size=$(csg_calc ${scale_factor} "*" $bin_size)
  table_begin=$(csg_calc ${scale_factor} "*" ${table_begin})
  table_end=$(csg_calc ${scale_factor} "*" ${table_end})
  scale="$(critical mktemp ${trunc}.pot.scale.XXXXX)"
  do_external table linearop --on-x "${input}" "${scale}" "${scale_factor}" "0"
else
  scale="${input}"
fi

#keep the grid for now, so that extrapolate can calculate the right mean
smooth="$(critical mktemp ${trunc}.pot.smooth.XXXXX)"
critical csg_resample --in ${scale} --out "$smooth" --grid "${table_begin}:${step}:${table_end}"

extrapol="$(critical mktemp ${trunc}.pot.extrapol.XXXXX)"
lfct="$(csg_get_interaction_property --allow-empty inverse.$sim_prog.table_left_extrapolation)"
rfct="$(csg_get_interaction_property --allow-empty inverse.$sim_prog.table_right_extrapolation)"
do_external potential extrapolate ${clean:+--clean} ${lfct:+--lfct ${lfct}} ${rfct:+--rfct ${rfct}} --type "$bondtype" "${smooth}" "${extrapol}"

interpol="$(critical mktemp ${trunc}.pot.interpol.XXXXX)"
deriv="$(critical mktemp ${trunc}.pot.deriv.XXXXX)"
critical csg_resample --in "${extrapol}" --out "$interpol" --grid "${table_begin}:${bin_size}:${table_end}" --der "${deriv}" --comment "$comment"

if [[ $do_shift = "yes" ]]; then
  tshift="$(critical mktemp ${trunc}.pot.shift.XXXXX)"
  do_external potential shift --type "$bondtype" "${interpol}" "${tshift}"
else
  tshift="$interpol"
fi

do_external convert_potential tab --header "${sim_prog}" --type "${bondtype}" "${tshift}" "${deriv}" "${output}"
if [[ $clean ]]; then
  rm -f "${smooth}" "${interpol}" "${extrapol}" "${tshift}"
  [[ ${input} != ${scale} ]] && rm -f "${scale}"
fi
