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
This script is a wrapper to convert a potential to gromacs

Usage: ${0##*/} [options] input output

Allowed options:
    --help       show this help
    --clean      remove all intermediate temp files
    --no-r2d     do not converts rad to degree (scale x axis with 180/3.1415)
                 for angle and dihedral
                 Note: VOTCA calcs in rad, but gromacs in degree
    --no-shift   do not shift the potential
    --step XXX   use XXX as step for the interaction
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
   --r2d) #default now
    shift ;;
   --no-r2d)
    r2d=1
    shift ;;
   --clean)
    clean="yes"
    shift ;;
   --no-shift)
    do_shift="no"
    shift ;;
   --step)
    step="$2"
    shift 2;;
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

#special if calling from csg_call
xvgtype="$(csg_get_interaction_property bondtype)"
[[ $xvgtype = "C6" || $xvgtype = "C12" || $xvgtype = "CB" ]] && tabtype="non-bonded" || tabtype="$xvgtype"

zero=0
if [[ $tabtype = "non-bonded" ]]; then
  tablend="$(csg_get_property --allow-empty cg.inverse.gromacs.table_end)"
  mdp="$(csg_get_property cg.inverse.gromacs.mdp)"
  if [[ -f ${mdp} ]]; then
    echo "Found setting file '$mdp' now trying to check options in there"
    rlist=$(get_simulation_setting rlist)
    tabext=$(get_simulation_setting table-extension)
    # if we have all 3 numbers do this checks
    tabl=$(csg_calc "$rlist" + "$tabext")
    [[ -n $tablend  ]] &&  csg_calc "$tablend" "<" "$tabl" && \
      die "${0##*/}: Error table is shorter then what mdp file ($mdp) needs, increase cg.inverse.gromacs.table_end in setting file.\nrlist ($rlist) + tabext ($tabext) > cg.inverse.gromacs.table_end ($tablend)"
    max="$(csg_get_interaction_property max)"
    rvdw="$(get_simulation_setting rvdw)"
    csg_calc "$max" ">" "$rvdw" && die "${0##*/}: rvdw ($rvdw) is smaller than max ($max)"
    [[ -z $tablend ]] && tablend=$(csg_calc "$rlist" + "$tabext")
  elif [[ -z $tablend ]]; then
    die "${0##*/}: cg.inverse.gromacs.table_end was not defined in xml seeting file"
  fi
elif [[ $tabtype = "bond" ]]; then
  tablend="$(csg_get_property cg.inverse.gromacs.table_end)"
elif [[ $tabtype = "angle" ]]; then
  tablend=180
elif [[ $tabtype = "dihedral" ]]; then
  zero="-180"
  tablend=180
else
  die "${0##*/}: Unknown interaction type $tabtype"
fi

[[ $step ]] || step=$(csg_get_interaction_property step)
gromacs_bins="$(csg_get_property cg.inverse.gromacs.table_bins)"
comment="$(get_table_comment $input)"

if [[ $tabtype = "angle" || $tabtype = "dihedral" ]] && [[ $r2d != 1 ]]; then
  scale="$(critical mktemp ${trunc}.pot.scale.XXXXX)"
  do_external table linearop --on-x "${input}" "${scale}" "$r2d" "0"
  step=$(csg_calc $r2d "*" $step)
else
  scale="${input}"
fi

#keep the grid for now, so that extrapolate can calculate the right mean
smooth="$(critical mktemp ${trunc}.pot.smooth.XXXXX)"
critical csg_resample --in ${scale} --out "$smooth" --grid "${zero}:${step}:${tablend}"

extrapol="$(critical mktemp ${trunc}.pot.extrapol.XXXXX)"
do_external potential extrapolate ${clean:+--clean} --type "$tabtype" "${smooth}" "${extrapol}"

interpol="$(critical mktemp ${trunc}.pot.interpol.XXXXX)"
critical csg_resample --in "${extrapol}" --out "$interpol" --grid "${zero}:${gromacs_bins}:${tablend}" --comment "$comment"

if [[ $do_shift = "yes" ]]; then
  tshift="$(critical mktemp ${trunc}.pot.shift.XXXXX)"
  do_external potential shift --type "$tabtype" "${interpol}" "${tshift}"
else
  tshift="$interpol"
fi

potmax="$(csg_get_property --allow-empty cg.inverse.gromacs.pot_max)"
do_external convert_potential xvg ${potmax:+--max} ${potmax} --type "${xvgtype}" "${tshift}" "${output}"
if [[ $clean ]]; then
  rm -f "${smooth}" "${interpol}" "${extrapol}" "${tshift}"
  [[ ${input} = ${scale} ]] || rm -f "${scale}"
fi
