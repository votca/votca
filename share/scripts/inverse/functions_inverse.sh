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
* log          = send a message to the logfile
* msg          = message to stdout and logfile
* die          = error message to stderr and logfile, 
                 and kills all csg process
* do_external  = get scriptname for sourcewrapper and run it
                 supports for_all
* for_all      = run a command for all non-bonded pairs
* logrun       = exec to log output
* run_or_exit  = logrun + die if error
* check_for    = checks if a binary exist in the path
* check_deps   = checks the dependencies of a script

Examples:
* log "Hi"
* msg "Hi"
* die "Error at line 99"
* do_external init gromacs NVT
* do_external init potential for_all bonded
* for_all bonded init_potential.sh 1 2 3
* logrun CMD 
* run_or_exit CMD

USES: \$CSGXMLFILE \$SOURCE_WRAPPER \$CSGLOG \$CSGRESTART csg_property

PROVIDES: log die msg csg_get_interaction_property csg_get_property csg_taillog do_external for_all is_done mark_done sed run_or_exit  printf

NEEDS:
EOF
exit 0
fi


log () {
  local warn
  if [ "$1" = "--no-warn" ]; then
    shift
    warn="no"
  else
    warn="yes"
  fi
  if [ -z "$LOG_REDIRECTED" ]; then
    if [ -n "$CSGLOG" ]; then
      echo -e "$*" >> $CSGLOG
    else
      echo -e "$*"
    fi
  else
    if [ "$warn" = "yes" ]; then
      echo -e "WARNING: Nested log call, when calling 'log $*'"
      echo -e "         log was redirected by '$LOG_REDIRECTED'"
      echo -e "         Try to avoid this, by removing one redirect, help: function_help"
    fi
    echo -e "$*" 
  fi
}
#echo a msg but log it too
msg() {
  [[ -z "$CSGLOG" ]] || log "$*"
  echo -e "$*"
}

unset -f die
die () {
  [[ -z "$CSGLOG" ]] || log --no-warn "$*"
  echo -e "$*" 1>&2
  log --no-warn "killing all processes...." 
  #send kill signal to all process within the process groups
  kill 0
  exit 1
}

#takes a task, find the according script and run it.
#first 2 argument are the task
do_external() {
  local script
  [[ -n "${SOURCE_WRAPPER}" ]] || die "do_external: SOURCE_WRAPPER is undefined"
  script="$($SOURCE_WRAPPER $1 $2)" || die "do_external: $SOURCE_WRAPPER $1 $2 failed" 
  shift 2
  #logrun do_external is a good combi to use
  log --no-warn "Running subscript '${script##*/} $*'"
  $script "$@" || die "do_external: $script $@ failed"
}

logrun(){
  local ret
  [[ -n "$1" ]] || die "logrun: missing argument"
  #--no-warn due to the fact that we get the warning anyway
  log --no-warn "logrun: run '$*'"
  if [ -z "$LOG_REDIRECTED" ]; then
    export LOG_REDIRECTED="logrun $*"
    if [ -n "$CSGLOG" ]; then 
      bash -c "$*" >> $CSGLOG 2>&1
      ret=$?
    else
      bash -c "$*" 2>&1
      ret=$?
    fi
    unset LOG_REDIRECTED
  else
    echo -e "WARNING: Nested log call, when calling 'logrun $*'"
    echo -e "         log was redirected by '$LOG_REDIRECTED'"
    echo -e "         Try to avoid this, by removing one redirect, help: function_help"
    bash -c "$*" 2>&1
    ret=$?
  fi
  return $ret
}

#useful subroutine check if a command was succesful AND log the output
run_or_exit() {
   logrun "$*" || die "run_or_exit: '$*' failed"
}

#do somefor all pairs, 1st argument is the type
for_all (){
  local bondtype csg_get
  local name interactions
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
  log "For all $bondtype"
  interactions="$(csg_property --file $CSGXMLFILE --short --path cg.${bondtype} --print name)" \
    || die "for_all: csg_property --file $CSGXMLFILE --short --path cg.${bondtype} --print name' failed"
  for name in $interactions; do
    log "for_all: run '$*'"
    #write variable defines in the front is better, that export
    #no need to run export -n afterwards
    bondtype="$bondtype" \
    bondname="$name" \
    bash -c "$*" || die "for_all: bash -c '$*' failed"   
  done
}

csg_taillog () {
  sync
  [[ -z "$CSGLOG" ]] || tail $* $CSGLOG
}

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
  cmd="csg_property --file $CSGXMLFILE --short --path cg.${bondtype} --filter \"name=$bondname\" --print $1"
  if ! ret="$($cmd 2>&1)"; then
    #csg_property failed and no default
    [[ -z "$2" ]] && die "csg_get_interaction_property: '$cmd' failed"
    ret="$2"
  fi
  [[ -z "$ret" ]] && [[ -n "$2" ]] && ret="$2"
  [[ "$allow_empty" = "no" ]] && [[ -z "$ret" ]] && \
    die "csg_get_interaction_property: Result of '$cmd' was empty"
  echo "$ret"
}

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
  cmd="csg_property --file $CSGXMLFILE --path '${1}' --short --print ."
  if ! ret="$($cmd 2>&1)"; then
    #csg_property failed and no default
    [[ -z "$2" ]] && die "csg_get_property: '$cmd' failed"
    ret="$2"
  fi
  [[ -z "$ret" ]] && [[ -n "$2" ]] && ret="$2"
  [[ "$allow_empty" = "no" ]] && [[ -z "$ret" ]] && \
    die "csg_get_property: Result of '$cmd' was empty"
  echo "$ret"
}

mark_done () {
  [[ -n "$1" ]] || die "mark_done: Missig argument"
  [[ -n "$CSGRESTART" ]] || die "mark_done: CSGRESTART is undefined"
  echo "$1 done" >> ${PWD}/$CSGRESTART 
}

is_done () {
  [[ -n "$1" ]] || die "is_done: Missig argument"
  [[ -n "$CSGRESTART" ]] || die "mark_done: CSGRESTART is undefined"
  [[ -f ${PWD}/${CSGRESTART} ]] || return 1
  [[ -n "$(sed -n "/^$1 done\$/p" ${PWD}/${CSGRESTART})" ]] && return 0
  return 1
}

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
    [[ -n "$(type -t $exe)" ]] || die "check_for: Could not find $exe needed ${file}" 
  done
}

check_deps () {
  [[ -n "$1" ]] || die "check_deps: Missig argument"
  local deps
  deps=$($1 --help | sed -n '/^USES:/p')
  [[ -z "${deps}" ]] && msg "check_for '$1' has no used block please add it" && return 0
  deps=$(echo "$deps" | sed 's/USES://')
  [[ -z "${deps}" ]] && return 0
  check_for "${1##*/}" $deps
}

int_check() {
  [[ -n "$2" ]] || die "int_check: Missig argument"
  [[ -z "${1//[0-9]}" ]] && return 0
  shift
  die "$*"
}

get_stepname() {
  local name
  [[ -n "$1" ]] || die "number_to_stepname: Missig argument"
  int_check "$1" "number_to_stepname needs a int as argument"
  name="$(printf step_%03i "$1")"
  [ -z "$name" ] && die "number_to_stepname: Could not get stepname"
  echo "$name"
}


#--------------------Exports-----------------
export -f die
export -f log 
export -f logrun 
export -f msg
export -f for_all
export -f run_or_exit
export -f do_external
export -f csg_get_property
export -f csg_get_interaction_property
export -f csg_taillog
export -f mark_done
export -f is_done
export -f check_for
export -f check_deps
export -f int_check
export -f get_stepname
