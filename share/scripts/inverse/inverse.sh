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
usage="Usage: ${0##*/} [OPTIONS] [setting_file.xml]"

show_help () {
  cat << eof
${0##*/}, version %version%



Start the script to run ibi, imc, etc.

$usage

Allowed options:
-h, --help                    show this help
-N, --do-iterations N         only do N iterations
    --wall-time SEK           Set wall clock time
    --options FILE            Specify the options xml file to use
    --clean                   clean out the PWD, dangerous

Examples:
* ${0##*/} cg.xml
* ${0##*/} -6 cg.xml

USES: csg_get_property date \$SOURCE_WRAPPER msg mkdir for_all do_external mark_done cp die is_done critical csg_get_interaction_property date \$CSGLOG date cp_from_main_dir get_current_step_dir get_last_step_dir get_main_dir get_stepname get_time rm update_stepnames

NEEDS: cg.inverse.method cg.inverse.program cg.inverse.iterations_max cg.inverse.filelist name cg.inverse.cleanlist
eof
}

#--help should always work so leave it here
if [ "$1" = "--help" ]; then
  show_help
  exit 0
fi

#do all start up checks option stuff
source "${0%/*}/start_framework.sh"  || exit 1

#defaults for options
do_iterations=""
wall_time=""

#unset stuff from enviorment
unset CSGXMLFILE CSGSCRIPTDIR CSGLOG

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
    int_check "$do_iterations" "inverse.sh: --do-iterations need a number as agrument"
    shift 2 ;;
   --wall-time)
    wall_time="$2"
    int_check "$wall_time" "inverse.sh: --wall-time need a number as agrument"
    start_time="$(date +%s)" || exit 1
    shift 2 ;;
   -[0-9]*)
    do_iterations=${1#-}
    shift ;;
   --clean)
    csg_ivnerse_clean
    exit $?;;
   --xmlfile)
    export CSGXMLFILE="$2"
    shift 2;;
   -h | --help)
    show_help
    exit 0;;
  *)
   die "Unknown option '$1'";;
 esac
done
### end parsing options

#old style maybe, new style set by --options
if [ -z "${CSGXMLFILE}" ]; then
  [ -n "$1" ] || die "Error: Missing xml file"
  export CSGXMLFILE="${1}"
  shift
fi

#make CSGXMLFILE a global path
[ "${CSGXMLFILE%/*}" = "${CSGXMLFILE}" ] && cpath="." || cpath="${CSGXMLFILE%/*}"
cpath="$(cd $cpath;pwd)"
export CSGXMLFILE="${cpath}/${CSGXMLFILE##*/}"
[ -f "$CSGXMLFILE" ] || die "${0##*/}: could not find ${CSGXMLFILE##*/} (long version: $CSGXMLFILE" 
unset cpath

#other stuff we need, which comes from xmlfile -> must be done here
#define $CSGRESTART
CSGRESTART="$(csg_get_property cg.inverse.restart_file)"
CSGRESTART="${CSGRESTART##*/}"
export CSGRESTART

#get csglog
CSGLOG="$(csg_get_property cg.inverse.log_file)"
CSGLOG="$PWD/${CSGLOG##*/}"
export CSGLOG
if [ -f "$CSGLOG" ]; then
  exec 3>&1 >> "$CSGLOG" 2>&1
  echo "\n\n#################################"
  echo "# Appending to existing logfile #"
  echo "#################################\n\n"
  echo "Sim started $(date)"
  msg "Appending to existing logfile ${CSGLOG##*/}"
else
  echo "For a more verbose log see: ${CSGLOG##*/}"
  #logfile is created in the next line
  exec 3>&1 >> "$CSGLOG" 2>&1
  echo "Sim started $(date)"
fi

method="$(csg_get_property cg.inverse.method)"
msg "We are doing Method: $method"

sim_prog="$(csg_get_property cg.inverse.program)"
echo "We are using Sim Program: $sim_prog"
source_function $sim_prog

iterations_max="$(csg_get_property cg.inverse.iterations_max)"
int_check "$do_iterations" "inverse.sh: cg.inverse.iterations_max needs to be a number"
echo "We are doing $iterations_max iterations (0=inf)."
convergence_check="$(csg_get_property cg.inverse.convergence_check "none")"
[ "$convergence_check" = "none" ] || echo "After every iteration we will do the following check: $convergence_check"

filelist="$(csg_get_property --allow-empty cg.inverse.filelist)"
[ -z "$filelist" ] || echo "We extra cp '$filelist' to every step to run the simulation"

cleanlist="$(csg_get_property --allow-empty cg.inverse.cleanlist)"
[ -z "$cleanlist" ] || echo "We extra clean '$cleanlist' after a step is done"

add_csg_scriptdir

critical $SOURCE_WRAPPER --status
critical $SOURCE_WRAPPER --check

#main script
[[ ! -f done ]] || { msg "Job is already done"; exit 0; }

update_stepnames 0
this_dir=$(get_current_step_dir --no-check)
if [ -d "$this_dir" ]; then
  msg "Skiping prepare"
  [[ -f $this_dir/done ]] || die "Incomplete step 0 (remove it if you don't know what to do)"
else
  msg ------------------------
  msg "Prepare (dir ${this_dir##*/})"
  msg ------------------------
  mkdir -p $this_dir || die "mkdir -p $this_dir failed"

  cd $this_dir || die "cd $this_dir failed"

  do_external prepare $method

  touch done
  msg "step 0 done"
  cd $(get_main_dir)
fi

begin=1
trunc=$(get_stepname --trunc)
for i in ${trunc}*; do
  [ -d "$i" ] || continue
  nr=${i#$trunc}
  if [ -n "$nr" ] && [ -z "${nr//[0-9]}" ]; then
    #convert to base 10, otherwise 008 is interpreted as octal
    nr=$((10#$nr))
    [ $nr -gt $begin ] && begin="$nr"
  fi
done
unset nr trunc
[ $begin -gt 1 ] && msg "Jumping in at iteration $begin"

avg_steptime=0
steps_done=0
[ $iterations_max -eq 0 ] && iterations=$begin || iterations=$iterations_max
for ((i=$begin;i<$iterations+1;i++)); do
  [ $iterations_max -eq 0 ] && ((iterations++))
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
      [[ -f "${this_dir}/${CSGRESTART}" ]] || die "No restart file found (remove this step if you don't know what to do - you will lose one step at max)"
    fi
  else
    echo "Step $i started at $(date)"
    mkdir -p $this_dir || die "mkdir -p $this_dir failed"
  fi

  cd $this_dir || die "cd $this_dir failed"
  is_done "stepdir" || mark_done "stepdir"

  if is_done "Initialize"; then
    msg "Initialization already done"
  else
    #get need files
    cp_from_main_dir $filelist

    #get files from last step, init sim_prog and ...
    do_external initstep $method

    mark_done "Initialize"
  fi

  if is_done "Simulation"; then
    msg "Simulation is already done"
  else
    msg "Simulation with $sim_prog"
    do_external run $sim_prog
    mark_done "Simulation"
  fi

  msg "Make update for $method"
  do_external update $method

  msg "Post update for $method"
  do_external post_update $method

  msg "Adding up potential for $method"
  do_external add_pot $method

  msg "Post add"
  do_external post add

  msg "Clean up"
  for cleanfile in ${cleanlist}; do
    rm -f $cleanfile
  done
  unset cleanfile

  step_time="$(( $(get_time) - $step_starttime ))"
  msg "\nstep $i done, needed $step_time secs"
  ((steps_done++))

  touch "done"

  if [ "$convergence_check" = "none" ]; then
    echo "No convergence check to be done"
  else
    msg "Doing convergence check: $convergence_check"
    if [ "$(do_external convergence_check "$convergence_check")" = "stop" ]; then
      msg "Iterations are converged, stopping"
      exit 0
    else
      msg "Iterations are not converged, going on"
    fi
  fi

  if [ -n "$wall_time" ]; then
    avg_steptime="$(( ( ( $steps_done-1 ) * $avg_steptime + $step_time ) / $steps_done + 1 ))"
    echo "New average steptime $avg_steptime"
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

touch "done"
echo "All done at $(date)"
exit 0

