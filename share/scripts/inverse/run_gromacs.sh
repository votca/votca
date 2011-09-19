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
This script runs a gromacs simulation or pre-simulation

Usage: ${0##*/} [--pre]

Used external packages: gromacs

EOF
   exit 0
fi

tpr="$(csg_get_property cg.inverse.gromacs.topol_out "topol.tpr")"

mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
[[ -f $mdp ]] || die "${0##*/}: gromacs mdp file '$mdp' not found"

conf="$(csg_get_property cg.inverse.gromacs.conf "conf.gro")"
[[ -f $conf ]] || die "${0##*/}: gromacs initial configuration file '$conf' not found"

confout="$(csg_get_property cg.inverse.gromacs.conf_out "confout.gro")"
mdrun_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.mdrun.opts)"

index="$(csg_get_property cg.inverse.gromacs.index "index.ndx")"
[[ -f $index ]] || die "${0##*/}: grompp index file '$index' not found"
topol="$(csg_get_property cg.inverse.gromacs.topol "topol.top")"
[[ -f $topol ]] || die "${0##*/}: grompp topol file '$topol' not found"

grompp_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"

grompp="$(csg_get_property cg.inverse.gromacs.grompp.bin "grompp")"
[ -n "$(type -p $grompp)" ] || die "${0##*/}: grompp binary '$grompp' not found"

ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
[[ $ext = "xtc" || $ext = "trr" ]] || die "${0##*/}: error trajectory type $ext is not supported"
if [[ $1 = "--pre" ]]; then
  : #in a presimulation usually do care about traj
elif  [[ $ext == "xtc" ]]; then
  [[ $(get_simulation_setting nstxtcout 0) -eq 0 ]] && die "${0##*/}: trajectory type (cg.inverse.gromacs.traj_type) is $ext, but nstxtcout is 0 in $mdp. Please check the setting again and remove the current step."
elif [[ $ext == "trr" ]]; then
  [[ $(get_simulation_setting nstxout 0) -eq 0 ]] && die "${0##*/}: trajectory type (cg.inverse.gromacs.traj_type) is $ext, but nstxout is 0 in $mdp. Please check the setting again and remove the current step."
fi

checkpoint="$(csg_get_property cg.inverse.gromacs.mdrun.checkpoint "state.cpt")"

#we want to do a presimulation and we are not currently performing it
if [[ $(csg_get_property cg.inverse.gromacs.pre_simulation "no") = "yes" && $1 != "--pre" ]]; then
  for i in tpr mdp conf confout index topol checkpoint; do
    f=${!i}
    [[ $f = */* ]] && die "${0##*/}: presimulation feature only work with local file (without /) in $i variable. Just try to copy $f to the maindir and add it to cg.inverse.filelist."
  done
  if ! is_done "Presimulation"; then #prepare/do presimulation
    [[ -d pre_simulation ]] || critical mkdir pre_simulation
    critical cp "${conf}" ./pre_simulation
    # a bit hacky but we have no better solution yet
    critical cp table*.xvg ./pre_simulation
    cp=0
    for i in mdp topol index; do
      f="$(csg_get_property "cg.inverse.gromacs.pre_simulation.$i" "${!i}")"
      [[ $f != ${!i} ]] && ((cp++))
      critical cp "${f}" "./pre_simulation/${!i}"
    done
    [[ $cp -eq 0 ]] && die "${0##*/}: mdp, topol and index of the presimulation are the same as for the main simulation, that does not make sense - at least one has to be different!"
    cd pre_simulation || die "${0##*/}: cd pre_simulation failed"
    msg "Doing pre simulation"
    do_external presimulation gromacs #easy to overwrite
    simulation_finish --no-traj || exit 0
    critical cp "${confout}" ../"${conf}"
    cd .. || die "${0##*/}: cd .. failed"
    mark_done "Presimulation" #after 'cd..' as restart file is here
  fi
  msg "Doing main simulation"
fi

critical $grompp -n "${index}" -f "${mdp}" -p "$topol" -o "$tpr" -c "${conf}" ${grompp_opts}
[ -f "$tpr" ] || die "${0##*/}: gromacs tpr file '$tpr' not found after runing grompp"

mdrun="$(csg_get_property cg.inverse.gromacs.mdrun.command "mdrun")"
#no check for mdrun, because mdrun_mpi could maybe exist only computenodes

if [ -n "$CSGENDING" ]; then
  #seconds left for the run
  wall_h=$(( $CSGENDING - $(get_time) ))
  #convert to hours
  wall_h=$(csg_calc $wall_h / 3600 )
  echo "${0##*/}: Setting $mdrun maxh option to $wall_h (hours)"
  mdrun_opts="-cpi $checkpoint -maxh $wall_h ${mdrun_opts}"
else
  echo "${0##*/}: No walltime defined, so no time limitation given to $mdrun"
fi

critical $mdrun -s "${tpr}" -c "${confout}" -o traj.trr -x traj.xtc ${mdrun_opts}
