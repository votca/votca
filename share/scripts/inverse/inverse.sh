#! /bin/bash
# 
# Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

#defaults
usage="Usage: ${0##*/} [OPTIONS] setting_file.xml"
do_iterations=""
clean="no"
wall_time=""

show_help () {
  cat << eof
${0##*/}, version %version%


Start the script to run ibm, imc, etc.

$usage

Allowed options:
-h, --help                    show this help
-N, --do-iterations N         only do N iterations
    --wall-time SEK           Set wall clock time
    --clean                   clean out the PWD, dangerous

Examples:
* ${0##*/} cg.xml
* ${0##*/} -6 cg.xml

USES: csg_get_property date \$SOURCE_WRAPPER msg mkdir for_all do_external mark_done cp die is_done log run_or_exit csg_get_interaction_property date \$CSGLOG date

NEEDS: cg.inverse.method cg.inverse.program cg.inverse.iterations_max cg.inverse.filelist name cg.inverse.cleanlist
eof
}

### begin parsing options
shopt -s extglob
while [ "${1#-}" != "$1" ]; do
 if [ "${1#--}" = "$1" ] && [ -n "${1:2}" ]; then
    #short opt with arguments here: fc
    if [ "${1#-[fc]}" != "${1}" ]; then
       set -- "${1:0:2}" "${1:2}" "${@:2}"
    else
       set -- "${1:0:2}" "-${1:2}" "${@:2}"
    fi
 fi
 case $1 in 
   --do-iterations)
    do_iterations="$2"
    shift 2 ;;
   --wall-time)
    wall_time="$2"
    start_time="$(date +%s)" || exit 1
    shift 2 ;;
   -[0-9]*)
    do_iterations=${1#-}
    shift ;;
   --clean)
    clean="yes"
    shift ;;
   -h | --help)
    show_help
    exit 0;;
  *)
   echo "Unknown option '$1'"
   exit 1;;
 esac
done

### end parsing options 

#do all start up checks option stuff
source "${0%.sh}_start.sh"  "$@" || exit 1
#shift away xml file
shift 1

#----------------End of pre checking--------------------------------
if [ -f "$CSGLOG" ]; then
  log "\n\n#################################"
  log "# Appending to existing logfile #"
  log "#################################\n\n"
  log "Sim started $(date)"
  echo "Appending to existing logfile ${CSGLOG##*/}"
else
  echo "For a more verbose log see: ${CSGLOG##*/}"
  #log is created in the next line
  echo "Sim started $(date)" > $CSGLOG || exit 1
fi

method="$(csg_get_property cg.inverse.method)" 
msg "We are doing Method: $method"
if [ "$method" = "imc" ]; then
  msg "####################################################"
  msg "# WARNING multicomponent imc is still experimental #"
  msg "####################################################"
fi

sim_prog="$(csg_get_property cg.inverse.program)"
log "We using Sim Program: $sim_prog"
source $($SOURCE_WRAPPER functions $sim_prog) || die "$SOURCE_WRAPPER functions $sim_prog failed" 

iterations="$(csg_get_property cg.inverse.iterations_max)" 
log "We are doing $iterations iterations."

filelist="$(csg_get_property --allow-empty cg.inverse.filelist)"
[ -z "$filelist" ] || log "We extra cp '$filelist' to every step to run the simulation"

cleanlist="$(csg_get_property --allow-empty cg.inverse.cleanlist)"
[ -z "$cleanlist" ] || log "We extra clean '$cleanlist' after a step is done"

run_or_exit $SOURCE_WRAPPER --status
run_or_exit $SOURCE_WRAPPER --check

#main script
[[ ! -f done ]] || { msg "Job is already done"; exit 0; }

update_stepnames 0
this_dir=$(get_current_step_dir --no-check)
if [ -d "$this_dir" ]; then
  msg "Skiping prepare"
  [[ -f $this_dir/done ]] || die "Incomplete step 0"
else
  msg ------------------------
  msg "Prepare (dir ${this_dir##*/})"
  msg ------------------------
  mkdir -p $this_dir || die "mkdir -p $this_dir failed"

  cd $this_dir || die "cd $this_dir failed"

  #copy+resample all rdf in $this_dir
  for_all non-bonded do_external resample calc

  do_external init $method

  #get confout.gro
  do_external init $sim_prog 

  for_all non-bonded cp '$(csg_get_interaction_property name).pot.new $(get_main_dir)' 
  touch done
  msg "step 0 done"
  cd $(get_main_dir)
fi

begin=1
trunc=$(get_stepname --trunc)
for i in ${trunc}*; do
  nr=${i#$trunc}
  [ -d "$i" ] && [ -n "$nr" ] && [ -z "${nr//[0-9]}" ] && [ $nr -gt $begin ] && begin="$nr"
done
unset nr trunc
[ $begin -gt 1 ] && msg "Jumping in at iteration $begin"

avg_steptime=0
steps_done=0
for ((i=$begin;i<$iterations+1;i++)); do
  step_starttime="$(get_time)"
  update_stepnames $i
  last_dir=$(get_last_step_dir)
  this_dir=$(get_current_step_dir --no-check)
  msg -------------------------------
  msg "Doing iteration $i (dir ${this_dir##*/})"
  msg -------------------------------
  if [ -d $this_dir ]; then
    if [ -f $this_dir/done ]; then
      msg "step $i is already done - skipping"
      continue
    else
      msg "Incomplete step $i"
      [[ -f "${this_dir}/${CSGRESTART}" ]] || die "No restart file found"
    fi
  else
    log "Step $i started at $(date)"
    mkdir -p $this_dir || die "mkdir -p $this_dir failed"
  fi

  cd $this_dir || die "cd $this_dir failed"

  if is_done "Initialize"; then
    msg "Initialization already done"
  else
    #copy+resample all rdf in this_dir 
    for_all non-bonded do_external resample calc

    #get need files
    cp_from_main_dir $filelist

    #get file from last step and so on
    do_external initstep $method

    #convert potential in format for sim_prog
    for_all non-bonded do_external convert_potential $sim_prog

    #Run simulation maybe change to Espresso or whatever
    do_external prepare $sim_prog 
    mark_done "Initialize"
  fi

  if is_done "Simulation"; then
    msg "Simulation is already done"
  else
    msg "Simulation runs"
    do_external run $sim_prog 
    mark_done "Simulation"
  fi

  msg "Make update $method" 
  do_external update $method

  msg "Post update"
  do_external post update 

  msg "Adding up potential"
  do_external add_pot $method

  msg "Post add"
  do_external post add

  #copy latest results
  for_all non-bonded 'cp $(csg_get_interaction_property name).pot.new $(get_main_dir)'

  touch done

  msg "Clean up"
  for cleanfile in ${cleanlist} ${CSGRESTART}; do
    logrun rm -f $cleanfile
  done
  unset cleanfile

  step_time="$(( $(get_time) - $step_starttime ))"
  msg "\nstep $i done, needed $step_time secs"
  ((steps_done++))

  if [ -n "$wall_time" ]; then
    avg_steptime="$(( ( ( $steps_done-1 ) * $avg_steptime + $step_time ) / $steps_done + 1 ))"
    log "New average steptime $avg_steptime"
    if [ $(( $(get_time) + $avg_steptime )) -gt $(( $wall_time + $start_time )) ]; then
      msg "We will not manage another step, stopping"
      exit 0
    else
      msg "We can go for another $(( ( ${start_time} + $wall_time - $(get_time) ) / $avg_steptime - 1 )) steps"
    fi
  fi

  if [ -n "$do_iterations" ]; then
    if [ $do_iterations -ge $steps_done ] ; then
      msg "Stopping at step $i, user requested to take some rest after this amount of iterations"
      exit 0
    else
      msg "Going on for another $(( $do_iterations - $steps_done )) steps"
    fi
  fi
  cd $(get_main_dir) || die "cd $(get_main_dir) failed"
done

touch done
log "All done at $(date)"
exit 0

