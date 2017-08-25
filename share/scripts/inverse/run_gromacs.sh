#! /bin/bash
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

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
This script runs a gromacs simulation or pre-simulation

Usage: ${0##*/} [--pre]

Used external packages: gromacs

EOF
   exit 0
fi

tpr="$(csg_get_property cg.inverse.gromacs.topol)"

mdp="$(csg_get_property cg.inverse.gromacs.mdp)"
[[ -f $mdp ]] || die "${0##*/}: gromacs mdp file '$mdp' not found (make sure it is in cg.inverse.filelist)"

conf="$(csg_get_property cg.inverse.gromacs.conf)"
[[ -f $conf ]] || die "${0##*/}: gromacs initial configuration file '$conf' not found (make sure it is in cg.inverse.filelist)"

confout="$(csg_get_property cg.inverse.gromacs.conf_out)"
mdrun_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.mdrun.opts)"

index="$(csg_get_property cg.inverse.gromacs.index)"
[[ -f $index ]] || die "${0##*/}: grompp index file '$index' not found (make sure it is in cg.inverse.filelist)"

topol_in="$(csg_get_property cg.inverse.gromacs.topol_in)"
[[ -f $topol_in ]] || die "${0##*/}: grompp text topol file '$topol_in' not found"

grompp_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"

grompp=( $(csg_get_property cg.inverse.gromacs.grompp.bin) )
[[ -n "$(type -p ${grompp[0]})" ]] || die "${0##*/}: grompp binary '${grompp[0]}' not found"

traj=$(csg_get_property cg.inverse.gromacs.traj)
if [[ $1 != "--pre" ]]; then
  #in a presimulation usually do care about traj and temperature
  if [[ $(csg_get_property cg.inverse.gromacs.runmode) = "none" ]]; then
  check_temp || die "${0##*/}: check of tempertures failed"
  fi
  if  [[ $traj == *.xtc ]]; then
    #XXX is returned if nstxout-compressed is not in mdp file
    nstxtcout=$(get_simulation_setting nstxout-compressed XXX)
    [[ ${nstxtcout} = XXX ]] && nstxtcout=$(get_simulation_setting nstxtcout 0)
    [[ ${nstxtcout} -eq 0 ]] && die "${0##*/}: trajectory type (cg.inverse.gromacs.traj) is '${traj##*.}', but nstxtcout is 0 in $mdp. Please check the setting again and remove the current step."
  elif [[ $traj == *.trr ]]; then
    [[ $(get_simulation_setting nstxout 0) -eq 0 ]] && die "${0##*/}: trajectory type (cg.inverse.gromacs.traj) is '${traj##*.}', but nstxout is 0 in $mdp. Please check the setting again and remove the current step."
  else
    die "${0##*/}: error trajectory type '${traj##*.}' (ending from '$traj') is not supported"
  fi
fi

checkpoint="$(csg_get_property cg.inverse.gromacs.mdrun.checkpoint)"

#we want to do a presimulation and we are not currently performing it
if [[ $(csg_get_property cg.inverse.gromacs.pre_simulation) = "yes" && $1 != "--pre" ]]; then
  for i in tpr mdp conf confout index topol_in checkpoint; do
    f=${!i}
    [[ $f = */* ]] && die "${0##*/}: presimulation feature only work with local file (without /) in $i variable. Just try to copy $f to the maindir and add it to cg.inverse.filelist."
  done
  if ! is_done "Presimulation"; then #prepare/do presimulation
    [[ -d pre_simulation ]] || critical mkdir pre_simulation
    critical cp "${conf}" ./pre_simulation
    # a bit hacky but we have no better solution yet
    critical cp table*.xvg ./pre_simulation
    cp=0
    for i in mdp topol_in index; do
      f="$(csg_get_property "cg.inverse.gromacs.pre_simulation.$i" "${!i}")" #filter me away
      [[ $f != ${!i} ]] && ((cp++))
      [[ -f ${f} ]] || die "${0##*/}: file '$f' not found (make sure it is in cg.inverse.filelist)"
      critical cp "${f}" "./pre_simulation/${!i}"
    done
    [[ $cp -eq 0 ]] && die "${0##*/}: mdp, topol_in and index of the presimulation are the same as for the main simulation, that does not make sense - at least one has to be different!"
    cd pre_simulation || die "${0##*/}: cd pre_simulation failed"
    msg "Doing pre simulation"
    do_external presimulation gromacs #easy to overwrite
    simulation_finish --no-traj || exit 0 #time out or crash - checkpoint_exist in inverse.sh will also look for a checkpoint in pre_simulation.
    critical cp "${confout}" ../"${conf}"
    cd .. || die "${0##*/}: cd .. failed"
    mark_done "Presimulation" #after 'cd..' as restart file is here
  fi
  msg "Doing main simulation"
fi

#support for older mdp file, cutoff-scheme = Verlet is default for >=gmx-5 now, but does not work with tabulated interactions
#XXX is returned if cutoff-scheme is not in mdp file
if [[ $(get_simulation_setting cutoff-scheme XXX) = XXX ]]; then
  echo "cutoff-scheme = Group" >> $mdp
  msg --color blue --to-stderr "Automatically added 'cutoff-scheme = Group' to $mdp, tabulated interactions only work with Group cutoff-scheme!"
fi

if [[ ${CSG_MDRUN_STEPS} ]]; then
  msg --color blue --to-stderr "Appending -nsteps ${CSG_MDRUN_STEPS} to mdrun options"
  mdrun_opts+=" -nsteps $CSG_MDRUN_STEPS"
fi

if [[ ${CSG_MDRUN_OPTS} ]]; then
  msg --color blue --to-stderr "Appending ${CSG_MDRUN_OPTS} to mdrun options"
  mdrun_opts+=" ${CSG_MDRUN_OPTS}"
fi

#support of REMD during gromacs simlation
if [[ $(csg_get_property cg.inverse.gromacs.runmode) = "replex" ]]; then
  PREFIX="$(csg_get_property cg.inverse.gromacs.runmode.replex.prefix)"
  read -a T <<<"$(csg_get_property cg.inverse.gromacs.runmode.replex.replicas)"
  j=0
  for i in "${T[@]}"
     do
       echo $i
       #TEMP has to be written as a place holer for the temperature in the *.mdp file
       sed "s/TEMP/$i/" grompp.mdp > grompp_$j.mdp
       critical ${grompp[@]} -n "${index}" -f "grompp_$j.mdp" -p "$topol_in" -o "$PREFIX"_"$j.tpr" -c "${conf}" ${grompp_opts} 2>&1 | gromacs_log "${grompp[@]} -n "${index}" -f "${mdp}" -p "$topol_in" -o "$tpr" -c "${conf}" ${grompp_opts}"
    #see can run grompp again as checksum of tpr does not appear in the checkpoint
       [[ -f "$PREFIX"_"$j.tpr" ]] || die "${0##*/}: gromacs tpr file '$tpr' not found after runing grompp"
       let j=$j+1
     done
  mdrun="$(csg_get_property cg.inverse.gromacs.runmode.replex.mdrun.command)"
fi

if [[ $(csg_get_property cg.inverse.gromacs.runmode) = "none" ]]; then
#see can run grompp again as checksum of tpr does not appear in the checkpoint
  critical ${grompp[@]} -n "${index}" -f "${mdp}" -p "$topol_in" -o "$tpr" -c "${conf}" ${grompp_opts} 2>&1 | gromacs_log "${grompp[@]} -n "${index}" -f "${mdp}" -p "$topol_in" -o "$tpr" -c "${conf}" ${grompp_opts}"
  [[ -f $tpr ]] || die "${0##*/}: gromacs tpr file '$tpr' not found after runing grompp"

  mdrun="$(csg_get_property cg.inverse.gromacs.mdrun.command)"
  #no check for mdrun, because mdrun_mpi could maybe exist only computenodes
fi

if [[ -n $CSGENDING ]]; then
  #seconds left for the run
  wall_h=$(( $CSGENDING - $(get_time) ))
  #convert to hours
  wall_h=$(csg_calc $wall_h / 3600 )
  echo "${0##*/}: Setting $mdrun maxh option to $wall_h (hours)"
  mdrun_opts="-cpi $checkpoint -maxh $wall_h ${mdrun_opts}"
else
  echo "${0##*/}: No walltime defined, so no time limitation given to $mdrun"
fi

#>gmx-5.1 has new handling of bonded tables
if [[ ${mdrun_opts} != *tableb* ]]; then
  tables=
  for i in table_[abd][0-9]*.xvg; do
    [[ -f $i ]] && tables+=" $i"
  done
  if [[ -n ${tables} ]]; then
    gmx_ver="$(critical ${grompp[@]} -h 2>&1)"
    shopt -s extglob
    [[ ${gmx_ver} = *"VERSION 5.1"+(.1|.2)* ]] && die "GROMACS 5.1 to 5.1.2 don't support tabulated bonds (http://redmine.gromacs.org/issues/1913), please update your GROMACS version" 
    msg --color blue --to-stderr "Automatically added '-tableb${tables} to mdrun options (add -tableb option to cg.inverse.gromacs.mdrun.opts yourself if this is wrong)"
    mdrun_opts+=" -tableb${tables}"
  fi
fi

if [[ $(csg_get_property cg.inverse.gromacs.runmode) = "none" ]]; then
critical $mdrun -s "${tpr}" -c "${confout}" -o "${traj%.*}".trr -x "${traj%.*}".xtc ${mdrun_opts} ${CSG_RUNTEST:+-v} 2>&1 | gromacs_log "$mdrun -s "${tpr}" -c "${confout}" -o "${traj%.*}".trr -x "${traj%.*}".xtc ${mdrun_opts}"
fi

#mdrun for REMD
if [[ $(csg_get_property cg.inverse.gromacs.runmode) = "replex" ]]; then
critical $mdrun ${mdrun_opts}
fi

[[ -z "$(sed -n '/[nN][aA][nN]/p' ${confout})" ]] || die "${0##*/}: There is a nan in '${confout}', this seems to be wrong."

