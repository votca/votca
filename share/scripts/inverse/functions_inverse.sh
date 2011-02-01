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

#-------------------defines----------------

if [ "$1" = "--help" ]; then
  cat <<EOF
${0##*/}, version %version%



We have defined some useful (?) functions:
* msg               = message to screen and logfile
* die               = error message to stderr and logfile,
                      kills all csg process
* do_external       = get scriptname for sourcewrapper and run it,
                      supports for_all
* for_all           = run a command for all non-bonded pairs
* critical          = run and die if error
* check_for         = checks if a binary exists in the path
* check_deps        = checks the dependencies of a script

Examples:
* echo "Hi"
* msg "Hi"
* die "Error at line 99"
* do_external init gromacs NVT
* do_external init potential for_all bonded
* for_all bonded init_potential.sh 1 2 3
* critical CMD

USES: \$CSGXMLFILE \$SOURCE_WRAPPER \$CSGLOG \$CSGRESTART csg_property printf cp date

PROVIDES: die msg csg_get_interaction_property csg_get_property do_external for_all is_done mark_done sed critical cat_external show_external check_for check_deps int_check get_stepname update_stepnames get_current_step_dir get_last_step_dir get_main_dir get_current_step_nr get_step_nr cp_from_to cp_from_main_dir cp_from_last_step get_time get_number_tasks

EOF
exit 0
fi

#echo a msg to the screen and send it to logfile too 
msg() {
  if [ "$1" = "--to-stderr" ]; then
    shift
    [ -n "$*" ] && echo -e "$*" >&2
  else
    [ -n "$*" ] && echo -e "$*"
  fi
  [ -n "${CSGLOG}" ] && [ -t 3 ] && echo -e "$*" >&3
}
export -f msg

unset -f die
die () {
  local pid pids c
  msg "#################################"
  msg "#################################"
  msg "# ERROR:                        #"
  msg "$*"
  msg "#################################"
  msg "#################################"
  [ -z "$CSGLOG" ] || msg "For details see $CSGLOG"
  if [ -n "${CSG_MASTER_PID}" ]; then
    #grabbing the pid group would be easier, but it would not work on AIX
    pid=$$
    pids="$$"
    c=0
    #find the parent of pid until we reach CSG_MASTER_PID
    until [ ${CSG_MASTER_PID} -eq $pid ]; do
      #get the parent pid using BSD style due to AIX
      pid=$(ps -o ppid= -p $pid 2>/dev/null)
      #store them in inverse order to kill parents before the child
      pids="$pid $pids"
      ((c++))
      #at max 100 iterations
      if [ $c -eq 10000 ]; then
        #failback to default, see comment below
        pids="0"
        break
      fi
    done
    if [ -n "${CSGLOG}" ]; then
      echo "die: (called from $$)  CSG_MASTER_PID is $CSG_MASTER_PID"
      echo "die: pids to kill: $pids"
    fi
    kill $pids
  else
    #send kill signal to all process within the process groups
    kill 0
  fi
  exit 1
}
export -f die

#takes a task, show content of the according script
cat_external() {
  local script
  [[ -n "${SOURCE_WRAPPER}" ]] || die "cat_external: SOURCE_WRAPPER is undefined"
  script="$($SOURCE_WRAPPER $1 $2)" || die "cat_external: $SOURCE_WRAPPER $1 $2 failed"
  cat "$script"
}
export -f cat_external

#takes a task, shows the according script
show_external() {
  local script
  [[ -n "${SOURCE_WRAPPER}" ]] || die "cat_external: SOURCE_WRAPPER is undefined"
  script="$($SOURCE_WRAPPER $1 $2)" || die "cat_external: $SOURCE_WRAPPER $1 $2 failed"
  echo "$script"
}
export -f show_external

#takes a task, find the according script and run it.
#first 2 argument are the task
do_external() {
  local script tags quiet="no"
  [ "$1" = "-q" ] && quiet="yes" && shift
  [[ -n "${SOURCE_WRAPPER}" ]] || die "do_external: SOURCE_WRAPPER is undefined"
  script="$($SOURCE_WRAPPER $1 $2)" || die "do_external: $SOURCE_WRAPPER $1 $2 failed"
  tags="$1 $2"
  shift 2
  [ "$quiet" = "no" ] && echo "Running subscript '${script##*/} $*'(from tags $tags)"
  $script "$@" 2>&1 || die "do_external: subscript $script $* (from tags $tags) failed"
}
export -f do_external

#useful subroutine check if a command was succesful AND log the output
critical() {
   "$@" || die "critical: '$*' failed"
}
export -f critical

#do somefor all pairs, 1st argument is the type
for_all (){
  local bondtype name names interactions i j
  if [ -z "$2" ]; then
    die "for_all need at least two arguments"
  fi
  bondtype="$1"
  shift
  #check that type is bonded or non-bonded
  if [ "$bondtype" != "non-bonded" ]; then
    die  "for_all: Argmuent 1 '$bondtype' is not non-bonded"
  fi
  [[ -n "$CSGXMLFILE" ]] || die "for_all: CSGXMLFILE is undefined"
  [[ -n "$(type -p csg_property)" ]] || die "for_all: Could not find csg_property"
  echo "For all $bondtype"
  interactions="$(csg_property --file $CSGXMLFILE --short --path cg.${bondtype} --print name)" \
    || die "for_all: csg_property --file $CSGXMLFILE --short --path cg.${bondtype} --print name' failed"
  names=( ${interactions} )
  for ((i=0;i<${#names[@]};i++)); do
    for ((j=i+1;j<${#names[@]};j++)); do
      [ "${names[$i]}" = "${names[$j]}" ] && die "for_all: the interaction name '${names[$i]}' appeared twice, this is not allowed"
    done
  done
  for name in $interactions; do
    echo "for_all: run '$*'"
    #we need to use bash -c here to allow things like $(csg_get_interaction_property xxx) in arguments
    #write variable defines in the front is better, that export
    #no need to run unset afterwards
    bondtype="$bondtype" \
    bondname="$name" \
    bash -c "$*" || die "for_all: bash -c '$*' failed for bondname '$name'"
  done
}
export -f for_all

#the save version of csg_get
csg_get_interaction_property () {
  local ret allow_empty cmd
  if [ "$1" = "--allow-empty" ]; then
    shift
    allow_empty="yes"
  else
    allow_empty="no"
  fi
  [[ -n "$1" ]] || die "csg_get_interaction_property: Missig argument"
  [[ -n "$CSGXMLFILE" ]] || die "csg_get_interaction_property: CSGXMLFILE is undefined"
  [[ -n "$bondtype" ]] || die "csg_get_interaction_property: bondtype is undefined"
  [[ -n "$bondname" ]] || die "csg_get_interaction_property: bondname is undefined"
  [[ -n "$(type -p csg_property)" ]] || die "csg_get_interaction_property: Could not find csg_property"
  cmd="csg_property --file $CSGXMLFILE --short --path cg.${bondtype} --filter name=$bondname --print $1"
  #the --filter option will make csg_property fail, don't stop if we have an default
  if ! ret="$($cmd 2>&1)"; then
    [ "$allow_empty" = "no" ] && [ -z "$2" ] && \
      die "csg_get_interaction_property:\n'$cmd'\nfailed geting '$1' with error msg:\n $ret\n and no default for $1"
    #ret has error message
    ret=""
  fi
  [ "$allow_empty" = "no" ] && [ -z "$ret" ] && [ -n "$2" ] && ret="$2"
  [[ "$allow_empty" = "no" ]] && [[ -z "$ret" ]] && \
    die "csg_get_interaction_property: Could not get '$1'\nResult of '$cmd' was empty"
  echo "$ret"
}
export -f csg_get_interaction_property

#get a property from xml
csg_get_property () {
  local ret allow_empty cmd
  if [ "$1" = "--allow-empty" ]; then
    shift
    allow_empty="yes"
  else
    allow_empty="no"
  fi
  [[ -n "$1" ]] || die "csg_get_property: Missig argument"
  [[ -n "$CSGXMLFILE" ]] || die "csg_get_property: CSGXMLFILE is undefined"
  [[ -n "$(type -p csg_property)" ]] || die "csg_get_property: Could not find csg_property"
  cmd="csg_property --file $CSGXMLFILE --path ${1} --short --print ."
  #csg_property only fails if xml file is bad otherwise result is empty
  ret="$(critical $cmd)"
  [[ -z "$ret" ]] && [[ -n "$2" ]] && ret="$2"
  [[ "$allow_empty" = "no" ]] && [[ -z "$ret" ]] && \
    die "csg_get_property: Could not get '$1'\nResult of '$cmd' was empty"
  echo "$ret"
}
export -f csg_get_property

mark_done () {
  [[ -n "$1" ]] || die "mark_done: Missig argument"
  [[ -n "$CSGRESTART" ]] || die "mark_done: CSGRESTART is undefined"
  is_done "$1" || echo "$1 done" >> ${PWD}/$CSGRESTART
}
export -f mark_done

is_done () {
  [[ -n "$1" ]] || die "is_done: Missig argument"
  [[ -n "$CSGRESTART" ]] || die "mark_done: CSGRESTART is undefined"
  [[ -f ${PWD}/${CSGRESTART} ]] || return 1
  [[ -n "$(sed -n "/^$1 done\$/p" ${PWD}/${CSGRESTART})" ]] && return 0
  return 1
}
export -f is_done

check_for () {
  [[ -n "$2" ]] || die "check_for: Missig arguments"
  file="$1"
  shift
  local exe
  for exe in $@; do
    if [ -z "${exe##\$*}" ]; then
      exe=${exe#\$}
      [[ -n "${!exe}" ]] || die "check_for: '${exe}' is undefined in ${file}"
      continue
    fi
    [[ -n "$(type -t $exe)" ]] || die "check_for: Could not find $exe needed by ${file}"
  done
}
export -f check_for

check_deps () {
  [[ -n "$1" ]] || die "check_deps: Missig argument"
  local deps
  deps="$(critical $1 --help)"
  deps="$(echo "$deps" | sed -n '/^USES:/p')" || die "check_deps: sed failed"
  [[ -z "${deps}" ]] && msg "check_for '$1' has no used block please add it" && return 0
  deps=$(echo "$deps" | sed 's/USES://')
  [[ -z "${deps}" ]] && return 0
  check_for "${1##*/}" $deps
}
export -f check_deps

int_check() {
  [[ -n "$2" ]] || die "int_check: Missig argument"
  [[ -z "${1//[0-9]}" ]] && return 0
  shift
  die "$*"
}
export -f int_check

get_stepname() {
  local name
  [[ -n "$1" ]] || die "get_stepname: Missig argument"
  if [ "$1" = "--trunc" ]; then
    echo "step_"
    return 0
  fi
  int_check "${1#-}" "get_stepname: needs a int as argument, but was $1"
  name="$(printf step_%03i "$1")"
  [ -z "$name" ] && die "get_stepname: Could not get stepname"
  echo "$name"
}
export -f get_stepname

update_stepnames(){
  local thisstep laststep nr
  [[ -n "$1" ]] || die "update_stepnames: Missig argument"
  nr="$1"
  int_check "$nr" "update_stepnames: needs a int as argument"
  [ -z "$CSG_MAINDIR" ] && die "update_stepnames: CSG_MAINDIR is empty"
  [ -d "$CSG_MAINDIR" ] || die "update_stepnames: $CSG_MAINDIR is not dir"
  thisstep="$(get_stepname $nr)"
  laststep="$(get_stepname $((nr-1)) )"
  export CSG_THISSTEP="$CSG_MAINDIR/$thisstep"
  export CSG_LASTSTEP="$CSG_MAINDIR/$laststep"
}
export -f update_stepnames

get_current_step_dir() {
  [ -z "$CSG_THISSTEP" ] && die "get_current_step_dir: CSG_THISSTEP is empty"
  if [ "$1" = "--no-check" ]; then
    :
  else
    [ -d "$CSG_THISSTEP" ] || die "get_last_step_dir: $CSG_THISSTEP is not dir"
  fi
  echo "$CSG_THISSTEP"

}
export -f get_current_step_dir

get_last_step_dir() {
  [ -z "$CSG_LASTSTEP" ] && die "get_last_step_dir: CSG_LASTSTEP is empty"
  [ -d "$CSG_LASTSTEP" ] || die "get_last_step_dir: $CSG_LASTSTEP is not dir"
  echo "$CSG_LASTSTEP"
}
export -f get_last_step_dir

get_main_dir() {
  [ -z "$CSG_MAINDIR" ] && die "get_main_dir: CSG_MAINDIR is empty"
  [ -d "$CSG_MAINDIR" ] || die "update_stepnames: $CSG_MAINDIR is not dir"
  echo "$CSG_MAINDIR"
}
export -f get_main_dir

get_current_step_nr() {
  local name nr
  name=$(get_current_step_dir)
  nr=$(get_step_nr $name)
  echo "$nr"
}
export -f get_current_step_nr

get_step_nr() {
  local nr trunc
  trunc=$(get_stepname --trunc)
  [[ -n "$1" ]] || die "get_step_nr: Missig argument"
  nr=${1##*/}
  nr=${nr#$trunc}
  #convert to base 10 and cut leading zeros
  nr=$((10#$nr))
  int_check "$nr" "get_step_nr: Could not fetch step nr"
  echo "$nr"
}
export -f get_step_nr

cp_from_to() {
  local i to from where
  if [ "$1" = "--from" ]; then
    from="$2"
    shift 2
  else
    die "cp_form_to: first argument should be --from DIR"
  fi
  if [ "$1" = "--where" ]; then
    where="$2"
    shift 2
  else
    where="."
  fi
  if [ "$1" = "--no-check" ]; then
    shift
  else
    [ -d "$where" ] || die "cp_from_to: $where is not a dir"
    [ -d "$from" ] || die "cp_from_to: $from is not a dir"
  fi
  [ -z "$1" ] && die "cp_from_to: Missing argument"
  for i in $@; do
    #no glob pattern in $i or could not be expanded
    if [ "$from/$i" = "$(echo $from/$i)" ]; then
      [ -e "$from/$i" ] || die "cp_from_to: could not find '$from/$i'"
    fi
    cp -r $from/$i "$where" 2>&1 || die "cp_from_to: cp -r '$from/$i' '$where' failed"
  done
}
export -f cp_from_to

cp_from_main_dir() {
  echo "cp_from_main_dir: '$@'"
  critical cp_from_to --from $(get_main_dir) "$@"
}
export -f cp_from_main_dir

cp_from_last_step() {
  echo "cp_from_last_step: '$@'"
  critical cp_from_to --from $(get_last_step_dir) "$@"
}
export -f cp_from_last_step

get_time() {
  date +%s || die "get_time:  time +%s failed"
}
export -f get_time

get_number_tasks() {
  local tasks
  [ -n "$(csg_get_property --allow-empty cg.inverse.mpi.tasks)" ] && \
    msg --to-stderr "get_number_tasks: the xml option cg.inverse.mpi.tasks has been renamed to cg.inverse.parallel.tasks\nPlease remove the obsolete cg.inverse.mpi block, it is not used anyway\n"
  tasks="$(csg_get_property cg.inverse.parallel.tasks 1)"
  [ "$tasks" = "auto" ] && tasks=0
  int_check "$tasks" "get_number_tasks: cg.inverse.parallel.tasks needs to be a number"
  #this only work for linux
  if [ $tasks -eq 0 ] && [ -r /proc/cpuinfo ]; then
    tasks=$(sed -n '/processor/p' /proc/cpuinfo | sed -n '$=')
    [[ -z "${tasks//[0-9]}" ]] || tasks=1
  fi
  [ $tasks -le 1 ] && tasks=1
  echo "$tasks"
}
export -f get_number_tasks

get_table_comment() {
  local version
  [[ -n "$(type -p csg_call)" ]] || die "get_defaults_comment: Could not find csg_version"
  version="$(csg_call --version)" || die "get_defaults_comment: csg_call --version failed"
  echo "Created on $(date) by $USER@$HOSTNAME"
  echo "called from $version" | sed "s/csg_call/${0##*/}/"
  echo "settings file: $CSGXMLFILE"
  echo "working directory: $PWD"
}
export -f get_table_comment

csg_ivnerse_clean() {
  echo -e "So, you want to clean?\n"
  echo "We will remove:"
  files="$(ls -d done ${CSGRESTART} ${CSGLOG##$PWD/} step_* *~ 2>/dev/null)"
  if [ -z "$files" ]; then
    echo "Nothing to clean"
  else
    echo $files
    echo -e "\nCTRL-C to stop it"
    for ((i=10;i>0;i--)); do
      echo -n "$i "
      sleep 1
    done
    rm -rf $files
    echo -e "\n\nDone, hope you are happy now"
  fi
}
export -f csg_ivnerse_clean

add_csg_scriptdir() {
  #make this work even if there is no xmlfile
  [ -n "${CSGXMLFILE}" ] && CSGSCRIPTDIR="$(csg_get_property --allow-empty cg.inverse.scriptdir)"
  #but mayit was define elsewhere
  if [ -n "$CSGSCRIPTDIR" ]; then
    #scriptdir maybe contains $PWD or something
    eval CSGSCRIPTDIR=$CSGSCRIPTDIR
    CSGSCRIPTDIR="$(cd $CSGSCRIPTDIR;pwd)"
    [[ -d "$CSGSCRIPTDIR" ]] || die "CSGSCRIPTDIR '$CSGSCRIPTDIR' is not a dir"
    export CSGSCRIPTDIR
    export PERL5LIB="$CSGSCRIPTDIR:$PERL5LIB"
  fi
}
export -f add_csg_scriptdir

source_function() {
  local function_file
  [ -n "$1" ] || die "source_function: Missig argument"
  [[ -n "${SOURCE_WRAPPER}" ]] || die "source_function: SOURCE_WRAPPER is undefined"
  function_file=$($SOURCE_WRAPPER functions $1) || die "source_function: $SOURCE_WRAPPER functions $1 failed"
  source ${function_file} || die "source_function: source ${function_file} failed"
}
export -f source_function
