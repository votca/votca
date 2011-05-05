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

get_from_mdp() { #gets a parameter (1st argument) from gromacs mdp file (2nd parameter)
  local res
  [[ -z $1 || -z $2 ]] && die "get_from_mdp: Missing argument (what file)"
  [[ -f $2 ]] || die "get_from_mdp: Could not read file '$2'"
  #1. strip comments
  #2. get important line
  #3. remove leading and tailing spaces
  res="$(sed -e '/^[[:space:]]*;/d' -e 's#;.*$##' "$2" | \
        sed -n -e "s#^[[:space:]]*$1[[:space:]]*=[[:space:]]*\(.*\)[[:space:]]*\$#\1#p" | \
	sed -e 's#^[[:space:]]*##' -e 's#[[:space:]]*$##')" || \
    die "get_from_mdp: sed failed"
  [[ -z "$res" ]] && [ -z "$3" ] && die "get_from_mdp: could not fetch $1 from $2, please add it"
  [ -n "$res" ] && echo "$res" || echo "$3"
}
export -f get_from_mdp

check_cutoff() { #compared current interactions cutoff vs rvdw, 
  local max rvdw res cutoff_check
  [[ -z $1 ]] && die "check_cutoff: Missing argument (mdp file)"
  cutoff_check=$(csg_get_property cg.inverse.gromacs.cutoff_check "yes")
  [ ${cutoff_check} = "no" ] && return 0
  max="$(csg_get_interaction_property max)"
  rvdw="$(get_from_mdp rvdw "$1")"
  csg_calc "$max" ">" "$rvdw" && die "Error in interaction '$bondname': rvdw ($rvdw) in $1 is smaller than max ($max)\n\
To ignore this check set cg.inverse.gromacs.cutoff_check to 'no'"
  return 0
}
export -f check_cutoff

check_temp() { #compares k_B T in xml with temp in mpd file
  local temp_check kbt temp res
  [[ -z $1 ]] && die "check_temp: Missing argument (mdp file)"
  temp_check=$(csg_get_property cg.inverse.gromacs.temp_check "yes")
  [ ${temp_check} = "no" ] && return 0
  #kbt in energy unit
  kbt="$(csg_get_property cg.inverse.kBT)"
  temp="$(get_from_mdp ref_t "$1")"
  #0.00831451 is k_b in gromacs untis see gmx manual chapter 2
  csg_calc "$kbt" "=" "0.00831451*$temp" || die "Error:  cg.inverse.kBT ($kbt) in xml seetings file differs from 0.00831451*ref_t ($temp) in $1\n\
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
