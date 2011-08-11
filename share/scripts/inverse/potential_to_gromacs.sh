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

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
This script is a wrapper to convert a potential to gromacs

Usage: ${0##*/} [--clean] input output

Allowed options:
    --help                    show this help
    --clean                   remove all intermediate temp files
EOF
  exit 0
fi

if [[ $1 = "--clean" ]]; then
  clean="yes"
  shift
else
  clean="no"
fi

if [[ -n $1 ]]; then
  name="${1%%.*}"
  input="$1"
  shift
else
  name=$(csg_get_interaction_property name)
  input="${name}.pot.cur"
fi
[[ -f $input ]] || die "${0##*/}: Could not find input file '$input'"

if [[ -n $1 ]]; then 
  output="$1"
  shift
else
  output="$(csg_get_interaction_property inverse.gromacs.table)"
fi

echo "Convert $input to $output"

zero=0
tabtype="$(csg_get_interaction_property bondtype)"
#do this with --allow-empty to avoid stoping if calling from csg_call
[[ $(csg_get_property --allow-empty cg.inverse.method) = "tf" ]] && tabtype="thermforce"

if [[ $tabtype = "non-bonded" || $tabtype = "C6" || $tabtype = "C12" ]]; then
  tablend="$(csg_get_property --allow-empty cg.inverse.gromacs.table_end)"
  mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
  if [[ -f ${mdp} ]]; then
    echo "Found setting file '$mdp' now trying to check options in there"
    rlist=$(get_simulation_setting rlist)
    tabext=$(get_simulation_setting table-extension)
    # if we have all 3 numbers do this checks
    tabl=$(csg_calc "$rlist" + "$tabext")
    [[ -n $tablend  ]] &&  csg_calc "$tablend" "<" "$tabl" && \
      die "${0##*/}: Error table is shorter then what mdp file ($mdp) needs, increase cg.inverse.gromacs.table_end in setting file.\nrlist ($rlist) + tabext ($tabext) > cg.inverse.gromacs.table_end ($tablend)"
    [[ -z $tablend ]] && tablend=$(csg_calc "$rlist" + "$tabext")
  elif [[ -z $tablend ]]; then
    die "${0##*/}: cg.inverse.gromacs.table_end was not defined in xml seeting file"
  fi
elif [[ $tabtype = "bonded" || $tabtype = "thermforce" ]]; then
  tablend="$(csg_get_property cg.inverse.gromacs.table_end)"
elif [[ $tabtype = "angle" ]]; then
  tablend=180
elif [[ $tabtype = "dihedral" ]]; then
  zero="-180"
  tablend=180
fi

gromacs_bins="$(csg_get_property cg.inverse.gromacs.table_bins)"
comment="$(get_table_comment $input)"

smooth="$(critical mktemp ${name}.pot.smooth.XXXXX)"
critical csg_resample --in ${input} --out "$smooth" --grid "${zero}:${gromacs_bins}:${tablend}" --comment "$comment"
extrapol="$(critical mktemp ${name}.pot.extrapol.XXXXX)"

tshift="$(critical mktemp ${name}.pot.shift.XXXXX)"
if [[ $tabtype = "non-bonded" || $tabtype = "C6" || $tabtype = "C12" ]]; then
  extrapol1="$(critical mktemp ${name}.pot.extrapol2.XXXXX)"
  do_external table extrapolate --function exponential --avgpoints 5 --region left "${smooth}" "${extrapol1}"
  do_external table extrapolate --function constant --avgpoints 1 --region right "${extrapol1}" "${extrapol}"
  do_external pot shift_nonbonded "${extrapol}" "${tshift}"
elif [[ $tabtype = "thermforce" ]]; then
  do_external table extrapolate --function constant --avgpoints 5 --region leftright "${smooth}" "${extrapol}"
  do_external pot shift_bonded "${extrapol}" "${tshift}"
else
  do_external table extrapolate --function exponential --avgpoints 5 --region leftright "${smooth}" "${extrapol}"
  do_external pot shift_bonded "${extrapol}" "${tshift}"
fi

potmax="$(csg_get_property --allow-empty cg.inverse.gromacs.pot_max)"
[[ -n ${potmax} ]] && potmax="--max ${potmax}"

do_external convert_potential xvg ${potmax} --type "${tabtype}" "${tshift}" "${output}"
[[ $clean = "yes" ]] && rm -f "${smooth}" "${extrapol}" "${tshift}" "${extrapol1}"
