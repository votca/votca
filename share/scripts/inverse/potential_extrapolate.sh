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

show_help () {
  cat <<EOF
${0##*/}, version %version%
This script extrapolates a potential in the correct way depending on its type.

Usage: ${0##*/} [options] input output

Allowed options:
    --help                    show this help
    --clean                   remove all intermediate temp files
    --type TYPE               type of the potential
                              possible:${pot_types}
    --avg-point INT           number of average points
                              default: $avg_points

EOF
}

clean="no"
pot_type="$2"
pot_types=" non-bonded bonded thermforce "
avg_points=5

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
   --clean)
    clean="yes"
    shift ;;
   --type)
    pot_type="$2"
  --avg-point)
    avg_points=$2;
    is_int "$2" || die "${0##*/}: argument of --avg-point should be int"
    shift 2;;
   -h | --help)
    show_help
    exit 0;;
  *)
   die "Unknown option '$1'";;
 esac
done
### end parsing options

[[ -z $pot_type ]] && die "${0##*/}: please specify add potential type (--type options) from:${pot_types}"
[[ -n ${pot_types//* $pot_type *} ]] && die "${0##*/}: given potential type is not in list${pot_types}"

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"

input="$1"
[[ -f $input ]] || die "${0##*/}: Could not find input file '$input'"

output="$1"

echo "Extrapolate $input to $output"

extrapol="$(critical mktemp ${name}.pot.extrapol.XXXXX)"
if [[ $pot_type = "non-bonded"  ]]; then
  intermediate="$(critical mktemp ${input}.onlyleft.XXXXX)"
  do_external table outputate --function exponential --avgpoints $avg_points --region left "${input}" "${intermediate}"
  do_external table outputate --function constant --avgpoints 1 --region right "${intermediate}" "${output}"
elif [[ $pot_type = "thermforce" ]]; then
  do_external table outputate --function constant --avgpoints $avg_points --region leftright "${input}" "${output}"
elif [[ $pot_type = "bonded"  || $pot_type = "angle" || $pot_type = "dihedral" ]]; then
  do_external table outputate --function exponential --avgpoints $avg_points --region leftright "${input}" "${output}"
else
  [[ -n ${pot_types//* $pot_type *} ]] && die "${0##*/}: given potential type is not in list${pot_types}"
fi

if [[ $clean = "yes" ]]; then
  rm -f "${intermediate}"
fi
