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

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%

Useful functions for gromacs:
EOF
sed -n 's/^\(.*\)([)] {[^#]*#\(.*\)$/* \1    -- \2/p' ${0}

echo
echo Used external packages: gromacs
  exit 0
fi

#check performed when this file is sourced
gmxrc="$(csg_get_property --allow-empty cg.inverse.gromacs.gmxrc)"
if [ -n "${gmxrc}" ]; then
  [ -f "$gmxrc" ] || die "functions_gromacs: Could not find gmxrc from xml file '$gmxrc'"
  msg "sourcing gromacs setting file: $gmxrc"
  source $gmxrc || die "functions_gromacs: error when 'source $gmxrc'"
fi
unset gmxrc

get_simulation_setting() { #gets a parameter (1st argument) from gromacs mdp file (2nd parameter)
  local res
  if [[ $1 = "--file" ]]; then
    mdp="$2"
    shift 2
  else
    mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
  fi
  [[ -z $1 ]] && die "get_simulation_setting: Missing argument (property)"
  [[ -f $mdp ]] || die "get_simulation_setting: Could not read setting file '$mdp'"
  #1. strip comments
  #2. get important line
  #3. remove leading and tailing spaces
  res="$(sed -e '/^[[:space:]]*;/d' -e 's#;.*$##' "$mdp" | \
        sed -n -e "s#^[[:space:]]*$1[[:space:]]*=\(.*\)\$#\1#p" | \
	sed -e 's#^[[:space:]]*##' -e 's#[[:space:]]*$##')" || \
    die "get_simulation_setting: sed failed"
  [[ -z $res && -z $2 ]] && die "get_simulation_setting: could not fetch $1 from $mdp and no default given, please add it in there"
  [[ -n $res ]] && echo "$res" || echo "$2"
}
export -f get_simulation_setting

check_cutoff() { #compared current interactions cutoff vs rvdw, 
  local max rvdw
  [[ "$(csg_get_property cg.inverse.gromacs.cutoff_check "yes")" = "no" ]] && return 0
  max="$(csg_get_interaction_property max)"
  rvdw="$(get_simulation_setting rvdw)"
  csg_calc "$max" ">" "$rvdw" && die "Error in interaction '$bondname': rvdw ($rvdw) in $1 is smaller than max ($max)\n\
To ignore this check set cg.inverse.gromacs.cutoff_check to 'no'"
  return 0
}
export -f check_cutoff

check_temp() { #compares k_B T in xml with temp in mpd file
  local kbt kbt2 temp
  [[ "$(csg_get_property cg.inverse.gromacs.temp_check "yes")" = "no" ]] && return 0
  #kbt in energy unit
  kbt="$(csg_get_property cg.inverse.kBT)"
  temp="$(get_simulation_setting ref_t)"
  #0.00831451 is k_b in gromacs untis see gmx manual chapter 2
  kbt2=$(csg_calc "$temp" "*" 0.00831451)
  csg_calc "$kbt" "=" "$kbt2" || die "Error:  cg.inverse.kBT ($kbt) in xml seetings file differs from 0.00831451*ref_t ($temp) in $1\n\
To ignore this check set cg.inverse.gromacs.temp_check to 'no'"
  return 0
}
export -f check_temp

simulation_finish() { #checks if simulation is finished
  local ext traj confout
  ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
  traj="traj.${ext}"
  confout="$(csg_get_property cg.inverse.gromacs.conf_out "confout.gro")"
  [ -f "$traj" ] && [ -f "$confout" ] && return 0
  return 1
}
export -f simulation_finish

checkpoint_exist() { #check if a checkpoint exists
  local checkpoint
  checkpoint="$(csg_get_property cg.inverse.gromacs.mdrun.checkpoint "state.cpt")"
  [ -f "$checkpoint" ] && return 0
  return 1
}
export -f checkpoint_exist

calc_begin_time() { #return the max of dt*frames and eqtime
  local dt equi_time first_frame
  dt=$(get_simulation_setting dt)
  first_frame="$(csg_get_property cg.inverse.gromacs.first_frame 0)"
  equi_time="$(csg_get_property cg.inverse.gromacs.equi_time 0)"
  t1=$(csg_calc "$dt" "*" "$first_frame")
  csg_calc "$t1" '>' "$equi_time" && echo "$t1" || echo "$equi_time"
}
export -f calc_begin_time

calc_end_time() { #return dt * nsteps
  local dt steps
  dt=$(get_simulation_setting dt)
  steps=$(get_simulation_setting nsteps)
  csg_calc "$dt" "*" "$steps"
}
export -f calc_end_time
