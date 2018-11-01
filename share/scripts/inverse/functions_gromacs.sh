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

Useful functions for gromacs:
EOF
sed -n 's/^\(.*\)([)] {[^#]*#\(.*\)$/* \1    -- \2/p' ${0}

echo
echo Used external packages: gromacs
  exit 0
fi

#check performed when this file is sourced
gmxrc="$(csg_get_property --allow-empty cg.inverse.gromacs.gmxrc)"
if [[ -n ${gmxrc} ]]; then
  [[ -f $gmxrc ]] || die "${BASH_SOURCE[0]##*/}: Could not find gmxrc from xml file '$gmxrc'"
  msg "sourcing gromacs setting file: $gmxrc"
  source $gmxrc || die "${BASH_SOURCE[0]##*/}: error when 'source $gmxrc'"
fi
unset gmxrc

#>gromacs-2016.3 allows to disable frame count output, e.g. "Reading frame ..."
export GMX_TRAJECTORY_IO_VERBOSITY=0

get_simulation_setting() { #gets a parameter (1st argument) from gromacs mdp file (default 2nd parameter)
  local res
  if [[ $1 = "--file" ]]; then
    mdp="$2"
    shift 2
  else
    mdp="$(csg_get_property cg.inverse.gromacs.mdp)"
  fi
  [[ -z $1 ]] && die "${FUNCNAME[0]}: Missing argument (property)"
  [[ -f $mdp ]] || die "${FUNCNAME[0]}: Could not read setting file '$mdp'"
  #1. strip comments
  #2. get important line
  #3. remove leading and tailing spaces
  res="$(sed -e '/^[[:space:]]*;/d' -e 's#;.*$##' "$mdp" | \
        sed -n -e "s#^[[:space:]]*$1[[:space:]]*=\(.*\)\$#\1#p" | \
	sed -e 's#^[[:space:]]*##' -e 's#[[:space:]]*$##')" || \
    die "${FUNCNAME[0]}: sed failed"
  [[ -z $res && -z $2 ]] && die "${FUNCNAME[0]}: could not fetch $1 from $mdp and no default given, please add it in there"
  [[ -n $res ]] && echo "$res" || echo "$2"
}
export -f get_simulation_setting

check_temp() { #compares k_B T in xml with temp in mpd file
  local kbt kbt2 temp t
  [[ "$(csg_get_property cg.inverse.gromacs.temp_check)" = "no" ]] && return 0
  #kbt in energy unit
  kbt="$(csg_get_property cg.inverse.kBT)"
  temp="$(get_simulation_setting "ref[_-]t")"
  for t in $temp; do
    #0.00831451 is k_b in gromacs untis see gmx manual chapter 2
    kbt2=$(csg_calc "$t" "*" 0.00831451)
    csg_calc "$kbt" "=" "$kbt2" || die "${FUNCNAME[0]}:  cg.inverse.kBT ($kbt) in xml seetings file differs from 0.00831451*ref_t ($t from $temp)\n\
To disable this check set cg.inverse.gromacs.temp_check to 'no'"
  done
  return 0
}
export -f check_temp

simulation_finish() { #checks if simulation is finished
  local traj confout
  traj=$(csg_get_property cg.inverse.gromacs.traj)
  confout="$(csg_get_property cg.inverse.gromacs.conf_out)"
  [[ $1 = "--no-traj" ]] && [[ -f $confout ]] && return 0
  [[ -f $traj ]] && [[ -f $confout ]] && return 0
  return 1
}
export -f simulation_finish

checkpoint_exist() { #check if a checkpoint exists
  local checkpoint
  checkpoint="$(csg_get_property cg.inverse.gromacs.mdrun.checkpoint)"
  [[ -f $checkpoint ]] && return 0
  [[ $(csg_get_property cg.inverse.gromacs.pre_simulation) = "yes" && -f pre_simulation/$checkpoint ]] && return 0
  return 1
}
export -f checkpoint_exist

calc_begin_time() { #return the max of dt*frames and eqtime
  local dt equi_time first_frame
  dt=$(get_simulation_setting dt)
  first_frame="$(csg_get_property cg.inverse.gromacs.first_frame)"
  equi_time="$(csg_get_property cg.inverse.gromacs.equi_time)"
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

gromacs_log() { #redirect stdin to a separate gromacs log file, 1st argument can be the name of the command to echo if redirection takes place
  local log log2
  log="$(csg_get_property --allow-empty cg.inverse.gromacs.log)"
  if [[ -z $log ]]; then
    [[ ${CSG_RUNTEST} ]] && tee >(cat - >&4) || cat
    return $?
  fi
  log="${log##*/}"
  log2="$(csg_get_property cg.inverse.log_file)"
  log2="${log2##*/}"
  [[ $log = $log2 ]] && die "${FUNCNAME}: cg.inverse.gromacs.log is equal cg.inverse.log_file"
  [[ -n $* ]] && echo "Sending output of '$*' to $log (also look for errors there)" || echo "Sending stdin to $log (also look for errors there)"
  [[ ${CSG_RUNTEST} ]] && tee >(cat - >&4) || cat >> "$log"
  return $?
}
export -f gromacs_log
