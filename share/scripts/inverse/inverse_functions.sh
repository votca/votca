#! /bin/bash

#-------------------defines----------------

#takes a task, find the according script and run it.
#first 2 argument are the task or --direct script
do_external() {
  local script
  script="$($SOURCE_WRAPPER $1 $2)" || return 1
  shift 2
  if [ "$1" = "for_all" ]; then
    bondtype="$2"
    shift 2
    echo Running subscript ${script##*/} "$@" for bondtypes $bondtype > /dev/stderr
    for_all "$bondtype" "$script" "$@" || return 1
  else
    echo Running subscript ${script##*/} "$@" > /dev/stderr
    "$script" "$@" || return 1
  fi
  return 0
}

#useful subroutine check if a command was succesful AND log the output
#if first argrument is --log, 2nd argument is log name
run_or_exit() {
   local prog mylog
   [[ "$1" = "--log" ]] && { mylog="$2"; shift 2; }
   prog=$1
   shift
   [[ -n "$prog" ]] || { echo Error give one argument >/dev/stderr; return 1; }
   [[ -z "$mylog" ]] && mylog="log_${prog##*/}"
   echo "Running ${prog##*/} $* &> $mylog" 
   $prog $* &> $mylog
   [[ $? -eq 0 ]] || { echo Error at ${prog##*/}; return 1; }
}

#do somefor all pairs, 1st argument is the type
for_all (){
  local bondtype type1 type2
  local values csg script
  script="yes"
  if [ -z "$2" ]; then
    echo for_all need at least two arguments > /dev/stderr
    return 1
  fi
  bondtype="$1"
  shift
  values="bonded non-bonded"
  #check that type is bonded or non-bonded
  if [ -z "$bondtype" ] || [ -n "${values//*$bondtype*}" ]; then
    echo "Argmuent 1 '$bondtype' is not in $values" > /dev/stderr
    return 1
  fi
  csg="csg_property --file $CSGXMLFILE --short --path cg.${bondtype}"
  for pair in $($csg --print "name"); do
    type1=$($csg --filter name="$pair" --print "type1")
    type2=$($csg --filter name="$pair" --print "type2")
    #write variable defines in the front is better, that export
    #no need to run export -n afterwards
    #$1 ${@:2} is a trick to make clear the $1 is the command
    name=$pair \
    bondtype=$bondtype \
    type1=$type1 \
    type2=$type2 \
    bash -xc "$1 ${@:2}" || return 1
  done
  return 0
}

#gets simulation property from xml
get_sim_property () {
  if [ -z "$1" ]; then
    echo Missig arrgument for get_sim_property > /dev/stderr
    return 1
  fi
  csg_property --file $CSGXMLFILE --path cg.inverse.${1} --short --print .
  return 0
}

#--------------------Exports-----------------
export -f for_all
export -f run_or_exit

