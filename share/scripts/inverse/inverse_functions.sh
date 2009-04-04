#! /bin/bash

#-------------------defines----------------

#takes a task, find the according script and run it.
#first 2 argument are the task or --direct script
do_external() {
  local script
  script="$($SOURCE_WRAPPER $1 $2)" || exit 1
  shift 2
  if [ "$1" = "for_all" ]; then
    bondtype="$2"
    shift 2
    echo Running subscript ${script##*/} "$@" for bondtypes $bondtypes
    for_all "$bondtype" "$script" "$@" || exit 1
  else
    echo Running subscript ${script##*/} "$@"
    "$script" "$@" || exit 1
  fi
}

#useful subroutine check if a command was succesful AND log the output
#if first argrument is --log, 2nd argument is log name
run_or_exit() {
   local prog mylog
   [[ "$1" = "--log" ]] && { mylog="$2"; shift 2; }
   prog=$1
   shift
   [[ -n "$prog" ]] || { echo Error give one argument >/dev/stderr; exit 1; }
   [[ -z "$mylog" ]] && mylog="log_${prog##*/}"
   echo "Running ${prog##*/} $* &> $mylog" 
   $prog $* &> $mylog
   [[ $? -eq 0 ]] || { echo Error at ${prog##*/}; exit 1; }
}

#do somefor all pairs, 1st argument is the type
for_all (){
  local bondtype type1 type2
  local values csg script
  script="yes"
  if [ -z "$2" ]; then
    echo for_all need at least two arguments > /dev/stderr
    exit 1
  fi
  bondtype="$1"
  shift
  values="bonded non-bonded"
  #check that type is bonded or non-bonded
  if [ -z "$bondtype" ] || [ -n "${values//*$bondtype*}" ]; then
    echo "Argmuent 1 '$bondtype' is not in $values" > /dev/stderr
    exit 1
  fi
  csg="csg_property --file $CSGXMLFILE --short --path cg.${bondtype}"
  for pair in $($csg --print "name"); do
    type1=$($csg --filter name="$pair" --print "type1")
    type2=$($csg --filter name="$pair" --print "type2")
    bondtype=$bondtype type1=$type1 type2=$type2 bash -xc "$@" || exit 1
  done
}

#--------------------Exports-----------------
export -f for_all
export -f run_or_exit
export -f do_external

