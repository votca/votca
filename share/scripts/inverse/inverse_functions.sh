#! /bin/bash

#takes a task, find the according script and run it.
#first 2 argument are the task or --direct script
do_external() {
   local script
   script="$($SOURCE_WRAPPER $1 $2)" || exit 1
   shift 2
   echo Running subscript ${script##*/} "$@"
   $script "$@"
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
   echo "Running ${prog##*/} $* &> $mylogi" 
   $prog $* &> $mylog
   [[ $? -eq 0 ]] || { echo Error at ${prog##*/}; exit 1; }
}

#do somefor all pairs, 1st argument is the type
#if 1st argument --cmd run a command instead of a script
for_all (){
  local type type1 type2
  local values csg script
  script="yes"
  if [ "$1" = "--cmd" ]; then
    shift
    script="no"
  fi
  if [ -z "$2" ]; then
    echo for_all need at least two arguments > /dev/stderr
    exit 1
  fi
  type="$1"
  values="bonded non-bonded"
  echo $CSGXMLFILE
  #check that type is bonded or non-bonded
  if [ -z "$type" ] || [ -n "${values//*$type*}" ]; then
    echo "Argmuent 1 '$type' is not in $values" > /dev/stderr
    exit 1
  fi
  csg="csg_property --file $CSGXMLFILE --short --path cg.$type"
  for pair in $($csg --print "name"); do
    type1=$($csg --filter name="$pair" --print "type1")
    type2=$($csg --filter name="$pair" --print "type2")
    if [ "$script" = "yes" ]; then
      echo type=$type type1=$type1 type2=$type2 "$@" || exit 1
    else
      echo eval "$@" || exit 1
    fi
  done
}

