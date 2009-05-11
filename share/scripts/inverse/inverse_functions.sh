#! /bin/bash

#-------------------defines----------------

log () {
  [[ -z "$LOG_REDIRECTED" ]] ||  die "Nested log call, when calling 'log $*',LOG_REDIRECTED was '$LOG_REDIRECTED',remove one loging !!!"
  export LOG_REDIRECTED="log $*" 
  echo -e "$*" >> $CSGLOG
  unset LOG_REDIRECTED
}
#echo a msg but log it too
msg() {
  log "$*"
  echo "$*"
}

unset -f die
die () {
  #same as log "$*", but avoid infinit log-die loop, when nested
  echo -e "$*" >> $CSGLOG
  echo "$*" 1>&2
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
  [[ -z "$LOG_REDIRECTED" ]] ||  die "Nested log call, when calling 'logrun $*', LOG_REDIRECTED was '$LOG_REDIRECTED', remove one loging!!!"
  log "logrun: run '$*'"
  export LOG_REDIRECTED="logrun $*" 
  bash -c "$*" >> $CSGLOG 2>&1
  ret=$?
  unset LOG_REDIRECTED
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

#gets simulation property from xml
get_sim_property () {
  if [ -z "$1" ]; then
    die "get_sim_property: Missig arrgument for get_sim_property" 
  fi
  csg_property --file $CSGXMLFILE --path cg.inverse.${1} --short --print . \
    || die "get_sim_property: csg_property --file $CSGXMLFILE --path cg.inverse.${1} --print . failed "
}

#get a property from xml
get_property () {
  if [ -z "$1" ]; then
    die "get_property: Missig arrgument for get_sim_property" 
  fi
  csg_property --file $CSGXMLFILE --path cg.${1} --short --print . \
    || die "get_property: csg_property --file $CSGXMLFILE --path cg.${1} --print . failed "
}

#--------------------Exports-----------------
export -f die
export -f log 
export -f logrun 
export -f msg
export -f for_all
export -f run_or_exit
export -f do_external
export -f get_sim_property
export -f get_property

