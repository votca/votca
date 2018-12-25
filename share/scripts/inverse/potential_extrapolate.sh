#!/bin/bash
#
# Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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
This script extrapolates a potential in the correct way depending on its type.

Usage: ${0##*/} [options] input output

Allowed options:
    --help                    show this help
    --clean                   remove all intermediate temp files
    --type TYPE               type of the potential
                              possible: ${pot_types}
    --lfct FCT                type of the left extrapolation function
                              possible: ${fct_type}
                              default: exponential(non-bonded), linear (bonded)
    --rfct FCT                type of the right extrapolation function
                              possible: ${fct_types}
                              default: constant(non-bonded), linear (bonded)
    --avg-point INT           number of average points
                              default: $avg_points

EOF
}

clean="no"
lfct=
rfct=
fct_types="constant linear quadratic exponential sasha"
pot_type=
pot_types="non-bonded bond angle dihedral"
avg_points=3

### begin parsing options
shopt -s extglob
while [[ ${1} = --* ]]; do
 case $1 in
   --clean)
    clean="yes"
    shift ;;
   --type)
    pot_type="$2"
    shift 2;;
   --[rl]fct)
    is_part "${2}" "${fct_types}" || die "${0##*/}: given function type ($2) is not in the list of functions (${fct_types})"
    eval "${1#--}"="$2"
    shift 2;;
  --avg-point)
    avg_points="$2";
    is_int "$2" || die "${0##*/}: argument of --avg-point should be int"
    shift 2;;
  --help)
    show_help
    exit 0;;
  *)
   die "Unknown option '$1'";;
 esac
done
### end parsing options

[[ -z $pot_type ]] && die "${0##*/}: please specify add potential type (--type options) from: ${pot_types}"
is_part "${pot_type}" "${pot_types}" || die "${0##*/}: given potential type($pot_type) is not in the list of types (${pot_types})"

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"

input="$1"
[[ -f $input ]] || die "${0##*/}: Could not find input file '$input'"

output="$2"

echo "Extrapolate $input to $output"

intermediate="$(critical mktemp "${input}.onlyleft.XXXXX")"
if [[ $pot_type = "non-bonded"  ]]; then
  do_external table extrapolate --function ${lfct:-exponential} --avgpoints $avg_points --region left "${input}" "${intermediate}"
  do_external table extrapolate --function ${rfct:-constant} --avgpoints 1 --region right "${intermediate}" "${output}"
elif [[ $pot_type = "bond"  || $pot_type = "angle" || $pot_type = "dihedral" ]]; then
  do_external table extrapolate --function ${lfct:-linear} --avgpoints $avg_points --region left "${input}" "${intermediate}"
  do_external table extrapolate --function ${rfct:-linear} --avgpoints $avg_points --region right "${intermediate}" "${output}"
else
  die "${0##*/}: I don't know how to extraploate potential type '$pot_type', go and implement it!"
fi

if [[ $clean = "yes" ]]; then
  rm -f "${intermediate}"
fi
