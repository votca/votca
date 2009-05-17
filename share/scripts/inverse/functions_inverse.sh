#! /bin/bash

#-------------------defines----------------

function_help() {
  cat <<EOF
  We have defined some useful (?) functions:
  log         = send a message to the logfile
  msg         = message to stdout and logfile
  die         = error message to stderr and logfile, 
                and kills all csg process
  do_external = get scriptname for sourcewrapper and run it
                supports for_all
  for_all     = run at command for all
  logrun      = exec to log output
  run_or_exit = logrun + die if error


  Examples:
    log "Hi"
    msg "Hi"
    die "Error at line 99"
    do_external init gromacs NVT
    do_external init potential for_all bonded
    for_all bonded init_potential.sh 1 2 3
    logrun mdrun
    run_or_exit mdrun
EOF
}

log () {
  if [ -z "$LOG_REDIRECTED" ]; then
    echo -e "$*" >> $CSGLOG
  else
    echo -e "WARNING: Nested log call, when calling 'log $*'"
    echo -e "         log was redirected by '$LOG_REDIRECTED'"
    echo -e "         Try to avoid this, by removing one redirect, help: function_help"
    echo -e "$*" 
  fi
}
#echo a msg but log it too
msg() {
  log "$*"
  echo -e "$*"
}

unset -f die
die () {
  #same as log "$*", but avoid infinit log-die loop, when nested
  log "$*"
  echo -e "$*" 1>&2
  log "killing all processes...." 
  #send kill signal to all process within the process groups
  kill 0
  exit 1
}

#takes a task, find the according script and run it.
#first 2 argument are the task
do_external() {
  local script bondtype
  script="$($SOURCE_WRAPPER $1 $2)" || die "do_external: $SOURCE_WRAPPER $1 $2 failed" 
  shift 2
  if [ "$1" = "for_all" ]; then
    bondtype="$2"
    shift 2
    log "Running subscript '${script##*/} $*' for bondtypes $bondtype"
    for_all "$bondtype" "$script" "$@" 
  else
    log "Running subscript '${script##*/} $*'"
    $script "$@" || die "do_external: $script $@ failed"
  fi
}

logrun(){
  local ret
  [[ -n "$1" ]] || die "logrun: missing argument"
  log "logrun: run '$*'"
  if [ -z "$LOG_REDIRECTED" ]; then
    export LOG_REDIRECTED="logrun $*" 
    bash -c "$*" >> $CSGLOG 2>&1
    ret=$?
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
  local pair pairs
  if [ -z "$2" ]; then
    die "for_all need at least two arguments"
  fi
  bondtype="$1"
  shift
  #check that type is bonded or non-bonded
  if [ "$bondtype" != "non-bonded" ]; then
    die  "for_all: Argmuent 1 '$bondtype' is not non-bonded" 
  fi
  csg_get="csg_property --file $CSGXMLFILE --short --path cg.${bondtype} --print"
  log "For all $bondtype"
  pairs="$($csg_get name)" || die "for_all: $csg --print name failed"
  for pair in $pairs; do
    #write variable defines in the front is better, that export
    #no need to run export -n afterwards
    log "for_all: run '$*'"
    bondtype=$bondtype \
    csg_get="$csg_get" \
    bash -c "$*" || die "for_all: bash -c '$*' failed"   
  done
}

csg_taillog () {
  sync
  tail $* $CSGLOG
}

#the save version of csg_get
csg_get () {
  local ret allow_empty
  if [ "$1" = "--allow-empty" ]; then
    shift
    allow_empty="yes"
  else
    allow_empty="no"
  fi
  [[ -n "$csg_get" ]] || die "csg_get: csg_get variable was empty"
  [[ -n "$1" ]] || die "csg_get: Missing argument"
  ret="$($csg_get $1)" || die "csg_get: '$csg_get $1' failed"
  [[ -n "$2" ]] && [[ -z "$ret" ]] && ret="$2"
  [[ "$allow_empty" = "no" ]] && [[ -z "$ret" ]] && \
    die "csg_get: Result of '$csg_get $1' was empty"
  echo "$ret"
}

#gets simulation property from xml
csg_get_sim_property () {
  local ret allow_empty
  if [ "$1" = "--allow-empty" ]; then
    shift
    allow_empty="yes"
  else
    allow_empty="no"
  fi
  [[ -n "$1" ]] || die "csg_get_sim_property: Missig argument" 
  ret="$(csg_property --file $CSGXMLFILE --path cg.inverse.${1} --short --print .)" \
    || die "csg_get_sim_property: csg_property --file $CSGXMLFILE --path cg.inverse.${1} --print . failed"
  [[ -n "$2" ]] && [[ -z "$ret" ]] && ret="$2"
  [[ "$allow_empty" = "no" ]] && [[ -z "$ret" ]] && \
    die "csg_get: Result of '$csg_get $1' was empty"
  echo "$ret"
}

#get a property from xml
csg_get_property () {
  local ret allow_empty
  if [ "$1" = "--allow-empty" ]; then
    shift
    allow_empty="yes"
  else
    allow_empty="no"
  fi
  [[ -n "$1" ]] || die "csg_get_property: Missig argument" 
  ret="$(csg_property --file $CSGXMLFILE --path cg.${1} --short --print .)" \
    || die "csg_get_property: csg_property --file $CSGXMLFILE --path cg.${1} --print . failed"
  [[ -n "$2" ]] && [[ -z "$ret" ]] && ret="$2"
  [[ "$allow_empty" = "no" ]] && [[ -z "$ret" ]] && \
    die "csg_get: Result of '$csg_get $1' was empty"
  echo "$ret"
}

#--------------------Exports-----------------
export -f die
export -f log 
export -f logrun 
export -f msg
export -f for_all
export -f run_or_exit
export -f do_external
export -f csg_get_sim_property
export -f csg_get_property
export -f csg_taillog
export -f function_help

