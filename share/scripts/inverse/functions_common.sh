#! /bin/bash
#
# Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

if [[ $1 = "--help" ]]; then
  cat <<EOF
${0##*/}, version %version%

This file defines some commonly used functions:
EOF
sed -n 's/^\(.*\)([)] {[^#]*#\(.*\)$/* \1  -- \2/p' ${0}
echo
exit 0
fi

export BASH #need in CsgFunctions.pm

shopt -s extglob

msg() { #echos a msg on the screen and send it to the logfile if logging is enabled
  local color colors="blue cyan cyann green red purp"
  if [[ -z ${CSGNOCOLOR} ]]; then
    local blue="[34;01m"
    local cyan="[36;01m"
    local cyann="[36m"
    local green="[32;01m"
    local red="[31;01m"
    local purp="[35;01m"
    local off="[0m"
  else
    local blue cyan cyann green red purp off
  fi
  if [[ $1 = "--color" ]]; then
    [[ -z $2 ]] && die "${FUNCNAME[0]}: missing argument after --color"
    is_part "$2" "${colors}" || die "${FUNCNAME[0]}: Unknown color ($colors allowed)"
    color="${!2}"
    shift 2
  else
    off=""
  fi
  if [[ $1 = "--to-stderr" ]]; then
    shift
    [[ -z $* ]] && return
    if [[ -n ${CSGLOG} ]]; then
      echo -e "${color}$*${off}" >&4
      echo -e "$*" >&2
    else
      echo -e "${color}$*${off}" >&2
    fi
  else
    [[ -z $* ]] && return
    if [[ -n ${CSGLOG} ]]; then
      echo -e "${color}$*${off}" >&3
      echo -e "$*"
    else
      echo -e "${color}$*${off}"
    fi
  fi
}
export -f msg

show_callstack() { #show the current callstack
  local space line
  if [[ -n $CSG_CALLSTACK ]]; then
    echo "$CSG_CALLSTACK"
    space="$(echo "$CSG_CALLSTACK" | sed -n '$s/[^[:space:]].*$/    /p')"
  else
    space=""
  fi
  [[ $0 = *"bash" ]] || echo "${space}${0} - linenumber ${BASH_LINENO[ $(( ${#FUNCNAME[@]} -2 ))]}"
  for ((c=${#FUNCNAME[*]}-1;c>0;c--)); do
    [[ ${FUNCNAME[$c]} = main ]] && continue #main is useless as the info was printed 2 lines above
    space+="    "
    if [[ $0 = *csg_call || $0 = *inverse.sh ]]; then
      echo "${space}${FUNCNAME[$c]} - linenumber ${BASH_LINENO[ $(( $c - 1 )) ]} in ${BASH_SOURCE[$c]}"
    else
      echo "${space}${FUNCNAME[$c]} - linenumber ${BASH_LINENO[ $(( $c - 1 )) ]} (see 'csg_call --cat function ${FUNCNAME[$c]}')"
    fi
  done
  [[ $1 = "--extra" ]] || return 0
  shift
  for i in "$@"; do
    space+="    "
    echo "${space}${i} - linenumber ???"
  done
}
export -f show_callstack

unset -f die
die () { #make the iterative frame work stopp
  local pid pids c place
  #Output callstack to stderr in case die was executed in $( )
  echo -e "\nCallstack:" >&2
  show_callstack >&2
  [[ -z $CSGLOG ]] && place="Details can be found above" || place="For details see the logfile $CSGLOG"
  msg --color red --to-stderr "$(csg_banner "ERROR:" "$@" "$place")"
  if [[ -n ${CSG_MASTER_PID} ]]; then
    #grabbing the pid group would be easier, but it would not work on AIX
    pid="$$"
    pids="$$"
    c=0
    #find the parent of pid until we reach CSG_MASTER_PID
    until [[ ${CSG_MASTER_PID} -eq $pid ]]; do
      #get the parent pid using BSD style due to AIX
      pid=$(ps -o ppid= -p "$pid" 2>/dev/null)
      [[ -z $pid ]] && pids="0" && break
      #store them in inverse order to kill parents before the child
      pids="$pid $pids"
      ((c++))
      #at max 10000 iterations
      if [[ $c -eq 10000 ]]; then
        #failback to default, see comment below
        pids="0"
        break
      fi
    done
    if [[ -n ${CSGLOG} ]]; then
      echo "${FUNCNAME[0]}: (called from $$)  CSG_MASTER_PID is $CSG_MASTER_PID" >&2
      echo "${FUNCNAME[0]}: pids to kill: $pids" >&2
    fi
    kill $pids
  else
    #send kill signal to all process within the process groups
    kill 0
  fi
  exit 1
}
export -f die

cat_external() { #takes a two tags and shows content of the according script
  local script
  script="$(source_wrapper $1 $2)" || die "${FUNCNAME[0]}: source_wrapper $1 $2 failed"
  if [[ $1 = "function" ]]; then
    type $2 | sed '1d'
  else
    cat "${script/ *}"
  fi
}
export -f cat_external

do_external() { #takes two tags, find the according script and excute it
  local script tags quiet="no" ret
  [[ $1 = "-q" ]] && quiet="yes" && shift
  script="$(source_wrapper $1 $2)" || die "${FUNCNAME[0]}: source_wrapper $1 $2 failed"
  tags="$1 $2"
  [[ $1 != "function" && ! -x ${script/ *} ]] && die "${FUNCNAME[0]}: subscript '${script/ *}' (from tags $tags), is not executable! (Run chmod +x ${script/ *})"
  #print this message to stderr to allow $(do_external ) and do_external XX > 
  [[ $quiet = "no" ]] && echo "Running subscript '${script##*/}${3:+ }${@:3}' (from tags $tags) dir ${script%/*}" >&2
  # in debugmode we don't need to do anything special for $1 = function as set -x is already done
  if [[ -n $CSGDEBUG ]] && [[ -n "$(sed -n '1s@bash@XXX@p' "${script/ *}")" ]]; then
    CSG_CALLSTACK="$(show_callstack)" "${BASH}" -x $script "${@:3}"
  elif [[ -n $CSGDEBUG && -n "$(sed -n '1s@perl@XXX@p' "${script/ *}")" ]]; then
    local perl_debug="$(mktemp perl_debug.XXX)" ret
    PERLDB_OPTS="NonStop=1 AutoTrace=1 frame=2 LineInfo=$perl_debug" perl -dS $script "${@:3}"
    ret=$?
    cat "$perl_debug" 2>&1
    [[ $ret -eq 0 ]]
  elif [[ $1 != "function" && -n "$(sed -n '1s@perl@XXX@p' "${script/ *}")" ]]; then
    CSG_CALLSTACK="$(show_callstack --extra "${script/ *}")" $script "${@:3}"
  else
    CSG_CALLSTACK="$(show_callstack)" $script "${@:3}"
  fi || die "${FUNCNAME[0]}: subscript" "$script ${*:3}" "(from tags $tags) failed"
}
export -f do_external

critical() { #executes arguments as command and calls die if not succesful
  local quiet="no"
  [[ $1 = "-q" ]] && quiet="yes" && shift
  [[ -z $1 ]] && die "${FUNCNAME[0]}: missing argument"
  #print this message to stderr because $(critical something) is used very often
  [[ $quiet = "no" ]] && echo "Running critical command '$*'" >&2
   "$@" || die "${FUNCNAME[0]}: '$*' failed"
}
export -f critical

for_all (){ #do something for all interactions (1st argument)
  local bondtype ibondtype rbondtype bondtypes name interactions quiet="no"
  [[ $1 = "-q" ]] && quiet="yes" && shift
  [[ -z $1 || -z $2 ]] && die "${FUNCNAME[0]}: need at least two arguments"
  bondtypes="$1"
  shift
  interactions=( $(csg_get_interaction_property --all name) )
  min=( $(csg_get_interaction_property --all min) )
  [[ ${#min[@]} -ne ${#interactions[@]} ]] && die "${FUNCNAME[0]}: one interaction has no name or min"
  name=$(has_duplicate "${interactions[@]}") && die "${FUNCNAME[0]}: interaction name $name appears twice"
  for bondtype in $bondtypes; do
    #check that type is bonded or non-bonded
    is_part "$bondtype" "non-bonded bonded angle bond dihedral" || die  "for_all: Argument 1 needs to be non-bonded, bonded, angle, bond or dihedral"
    [[ $quiet = "no" ]] && echo "For all $bondtype" >&2
    #internal bondtype
    is_part "$bondtype" "angle bond dihedral bonded" && ibondtype="bonded" || ibondtype="non-bonded"
    interactions=( $(csg_get_property --allow-empty cg.$ibondtype.name) ) #filter me away
    for name in "${interactions[@]}"; do
      #check if interaction is actually angle, bond or dihedral
      if is_part "$bondtype" "angle bond dihedral"; then
	rbondtype=$(bondtype="$ibondtype" bondname="$name" csg_get_interaction_property bondtype)
	[[ $rbondtype = $bondtype ]] || continue
      fi
      #print this message to stderr to avoid problem with $(for_all something)
      [[ $quiet = no ]] && echo "for_all: run '$*' for interaction named '$name'" >&2
      #we need to use bash -c here to allow things like $(csg_get_interaction_property name) in arguments
      #write variable defines in the front is better, that export
      #no need to run unset afterwards
      bondtype="$ibondtype" \
      bondname="$name" \
      CSG_CALLSTACK="$(show_callstack)" \
      "${BASH}" -c "$*" || die "${FUNCNAME[0]}: ${BASH} -c '$*' failed for interaction named '$name'"
    done
  done
}
export -f for_all

csg_get_interaction_property () { #gets an interaction property from the xml file, should only be called from inside a for_all loop or with --all option
  local ret allow_empty="no" for_all="no" xmltype
  while [[ $1 = --* ]]; do
    case $1 in
      --allow-empty)
        allow_empty="yes";;
      --all)
        for_all="yes";;
      *)
	die "${FUNCNAME[0]}: Unknow option '$1'";;
    esac
    shift
  done
  [[ -n $1 ]] || die "${FUNCNAME[0]}: Missing argument"

  if [[ $for_all = "yes" ]]; then
    [[ $1 = "bondtype" ]] && die "${FUNCNAME[0]}: --all + bondtype not implemented yet"
    local t
    for t in non-bonded bonded; do
      ret+=" $(csg_get_property --allow-empty "cg.$t.$1")" #filter me away
    done
    ret="$(echo "$ret" | trim_all)"
    [[ -z $ret ]] && die "${FUNCNAME[0]}: Not a single interaction has a value for the property $1"
    echo "$ret"
    return 0
  fi

  #make these this case work even without name or type (called by csg_call)
  if [[ $1 = "name" ]]; then
    [[ -n $bondname ]] && echo "$bondname" && return 0
    die "${FUNCNAME[0]}: bondname is undefined (when calling from csg_call set it by --ia-name option)"
  fi
  if [[ $1 = "bondtype" ]]; then
    #bondtype is special -> dirty hack - removed whenever issue 13 is fixed
    [[ -z "$bondtype" ]] && die "${FUNCNAME[0]}: bondtype is undefined (when calling from csg_call set it by --ia-type option)"
    #for_all notation for any kind of bonded interaction, find the real type
    if [[ $bondtype = "bonded" ]]; then
      [[ -z ${bondname} ]] && die "${FUNCNAME[0]}: bondtype 'bonded' needs a bondname (when calling from csg_call set it by --ia-name option) or change type to angle, bond or dihedral"
      [[ -n "$(type -p csg_property)" ]] || die "${FUNCNAME[0]}: Could not find csg_property"
      mapping="$(csg_get_property --allow-empty cg.inverse.map)" #make error message more useful
      [[ -z ${mapping} ]] && die "${FUNCNAME[0]}: bondtype 'bonded' needs a mapping file (cg.inverse.map in xml) to determine the actual bond type (when calling from csg_call better use --ia-type bond, angle or dihedral)"
      local map names=() ret= ret2 dup
      for map in ${mapping}; do
        [[ -f "$(get_main_dir)/$map" ]] || die "${FUNCNAME[0]}: Mapping file '$map' for bonded interaction not found in maindir"
	names+=( $(critical -q csg_property --file "$(get_main_dir)/$map" --path cg_molecule.topology.cg_bonded.*.name --print . --short) )
	[[ -n ${names[@]} ]] && dup=$(has_duplicate "${names[@]}") && die "${FUNCNAME[0]}: cg_bonded name '$dup' appears twice in file(s) $mapping"
        ret2="$(critical -q csg_property --file "$(get_main_dir)/$map" --path cg_molecule.topology.cg_bonded.* --filter name="$bondname" --print . --with-path | trim_all)"
        ret2="$(echo "$ret2" | critical sed -n 's/.*cg_bonded\.\([^[:space:]]*\) .*/\1/p')"
	if [[ -n $ret2 ]]; then
	  [[ -n $ret ]] && die "${FUNCNAME[0]}: Found cg_bonded type for name '$bondname' twice"
	  ret="${ret2}"
	fi
      done
      [[ -z $ret ]] && die "${FUNCNAME[0]}: Could not find a bonded definition with name '$bondname' in the mapping file(s) '$mapping'. Make sure to use the same name in the settings file (or --ia-name when calling from csg_call) and the mapping file."
      echo "$ret"
    else
      echo "$bondtype"
    fi
    return 0
  fi

  [[ -n "$CSGXMLFILE" ]] || die "${FUNCNAME[0]}: CSGXMLFILE is undefined (when calling from csg_call set it by --options option)"
  [[ -n $bondtype ]] || die "${FUNCNAME[0]}: bondtype is undefined (when calling from csg_call set it by --ia-type option)"
  [[ -n $bondname ]] || die "${FUNCNAME[0]}: bondname is undefined (when calling from csg_call set it by --ia-name option)"

  #map bondtype back to tags in xml file (for csg_call)
  case "$bondtype" in
    "non-bonded")
      xmltype="non-bonded";;
    "bonded"|"bond"|"angle"|"dihedral")
      xmltype="bonded";;
    *)
      msg "Unknown bondtype '$bondtype' - assume non-bonded"
      xmltype="non-bonded";;
  esac

  [[ -n "$(type -p csg_property)" ]] || die "${FUNCNAME[0]}: Could not find csg_property"
  #the --filter/--path(!=.) option will make csg_property fail if $1 does not exist
  #so no critical here
  ret="$(csg_property --file $CSGXMLFILE --short --path cg.${xmltype} --filter name=$bondname --print $1 | trim_all)"
  #overwrite with function call value
  [[ -z $ret && -n $2 ]] && ret="$2"
  [[ -z $ret ]] && echo "${FUNCNAME[0]}: No value for '$1' found in $CSGXMLFILE, trying ${VOTCA_CSG_DEFAULTS}" >&2
  # if still empty fetch it from defaults file
  if [[ -z $ret && -f ${VOTCA_CSG_DEFAULTS} ]]; then
    ret="$(critical -q csg_property --file "${VOTCA_CSG_DEFAULTS}" --short --path cg.${xmltype}.$1 --print . | trim_all)"
    [[ $allow_empty = "yes" && -n "$res" ]] && msg "WARNING: '${FUNCNAME[0]} $1' was called with --allow-empty, but a default was found in '${VOTCA_CSG_DEFAULTS}'"
    #from time to time the default is only given in the non-bonded section
    [[ -z $ret ]] && ret="$(critical -q csg_property --file "${VOTCA_CSG_DEFAULTS}" --short --path cg.non-bonded.$1 --print . | trim_all)"
    [[ -n $ret ]] && echo "${FUNCNAME[0]}: value for '$1' from ${VOTCA_CSG_DEFAULTS}: $ret" >&2
  fi
  [[ $allow_empty = "no" && -z $ret ]] && die "${FUNCNAME[0]}: Could not get '$1' for interaction with name '$bondname' from ${CSGXMLFILE} and no default was found in ${VOTCA_CSG_DEFAULTS}"
  [[ -z $ret ]] && echo "${FUNCNAME[0]}: returning emtpy value for '$1'" >&2
  echo "${ret}"
}
export -f csg_get_interaction_property

csg_get_property () { #get an property from the xml file
  local ret allow_empty
  if [[ $1 = "--allow-empty" ]]; then
    shift
    allow_empty="yes"
  else
    allow_empty="no"
  fi
  [[ -n $1 ]] || die "${FUNCNAME[0]}: Missing argument"
  [[ -n "$CSGXMLFILE" ]] || die "${FUNCNAME[0]}: CSGXMLFILE is undefined (when calling from csg_call set it by --options option)"
  [[ -n "$(type -p csg_property)" ]] || die "${FUNCNAME[0]}: Could not find csg_property"
  #csg_property only fails if xml file is bad otherwise result is empty
  #leave the -q here to avoid flooding with messages
  ret="$(critical -q csg_property --file "$CSGXMLFILE" --path ${1} --short --print . | trim_all)"
  #overwrite with function call value
  [[ -z $ret && -n $2 ]] && ret="$2"
  [[ -z $ret ]] && echo "${FUNCNAME[0]}: No value for '$1' found in $CSGXMLFILE, trying ${VOTCA_CSG_DEFAULTS}" >&2
  #if still empty fetch it from defaults file
  if [[ -z $ret && -f ${VOTCA_CSG_DEFAULTS} ]]; then
    ret="$(critical -q csg_property --file "${VOTCA_CSG_DEFAULTS}" --path "${1}" --short --print . | trim_all)"
    [[ $allow_empty = "yes" && -n "$res" ]] && msg "WARNING: '${FUNCNAME[0]} $1' was called with --allow-empty, but a default was found in '${VOTCA_CSG_DEFAULTS}'"
    #avoid endless recursion
    [[ $1 = cg.inverse.program && -n $ret ]] && sim_prog="$ret" || \
      sim_prog="$(csg_get_property cg.inverse.program)" #no problem to call recursively as sim_prog has a default
    if [[ -z $ret ]] && [[ $1 = *${sim_prog}* ]]; then
      local path=${1/${sim_prog}/sim_prog}
      ret="$(critical -q csg_property --file "${VOTCA_CSG_DEFAULTS}" --path "${path}" --short --print . | trim_all)"
    fi
    [[ -n $ret ]] && echo "${FUNCNAME[0]}: value for '$1' from ${VOTCA_CSG_DEFAULTS}: $ret" >&2
    [[ $allow_empty = "yes" && -n "$res" ]] && msg "WARNING: '${FUNCNAME[0]} $1' was called with --allow-empty, but a default was found in '${VOTCA_CSG_DEFAULTS}'"
  fi
  [[ $allow_empty = "no" && -z $ret ]] && die "${FUNCNAME[0]}: Could not get '$1' from ${CSGXMLFILE} and no default was found in ${VOTCA_CSG_DEFAULTS}"
  [[ -z $ret ]] && echo "${FUNCNAME[0]}: returning emtpy value for '$1'" >&2
  echo "${ret}"
}
export -f csg_get_property

trim_all() { #make multiple lines into one and strip white space from beginning and the end, reads from stdin
  [[ -n "$(type -p tr)" ]] || die "${FUNCNAME[0]}: Could not find tr"
  tr '\n' ' ' | sed -e s'/^[[:space:]]*//' -e s'/[[:space:]]*$//' || die "${FUNCNAME[0]}: sed of argument $i failed"
}
export -f trim_all

mark_done () { #mark a task (1st argument) as done in the restart file
  local file
  [[ -n $1 ]] || die "${FUNCNAME[0]}: Missing argument"
  file="$(get_restart_file)"
  is_done "$1" || echo "$1 done" >> "${file}"
}
export -f mark_done

is_done () { #checks if something is already do in the restart file
  local file
  [[ -n $1 ]] || die "${FUNCNAME[0]}: Missing argument"
  file="$(get_restart_file)"
  [[ -f ${file} ]] || return 1
  [[ -n "$(sed -n "/^$1 done\$/p" ${file})" ]] && return 0
  return 1
}
export -f is_done

is_int() { #checks if all arguments are integers
  local i
  [[ -z $1 ]] && die "${FUNCNAME[0]}: Missing argument"
  for i in "$@"; do
    [[ -n $i && -z ${i//[0-9]} ]] || return 1
  done
  return 0
}
export -f is_int

to_int() { #convert all given numbers to int using awk's int function
  local i
  [[ -z $1 ]] && die "${FUNCNAME[0]}: Missing argument"
  for i in "$@"; do
    is_num "$i" || die "${FUNCNAME[0]}: $i is not a number"
    awk -v x="$i" 'BEGIN{ print ( int(x) ) }' || die "${FUNCNAME[0]}: awk failed"
  done
  return 0
}
export -f to_int

is_part() { #checks if 1st argument is part of the set given by other arguments
  [[ -z $1 || -z $2 ]] && die "${FUNCNAME[0]}: Missing argument"
  [[ " ${@:2} " = *" $1 "* ]]
}
export -f is_part

has_duplicate() { #check if one of the arguments is double
  local i j
  [[ -z $1 ]] && die "${FUNCNAME[0]}: Missing argument"
  for ((i=1;i<$#;i++)); do
    for ((j=i+1;j<=$#;j++)); do
      [[ ${!i} = ${!j} ]] && echo ${!i} && return 0
    done
  done
  return 1
}
export -f has_duplicate

remove_duplicate() { #remove duplicates list of arguments
  local i j out=() c
  [[ -z $1 ]] && die "${FUNCNAME[0]}: Missing argument"
  for ((i=1;i<=$#;i++)); do
    c=0
    for ((j=0;j<${#out[@]};j++)); do
      [[ ${!i} = ${out[j]} ]] && ((c++))
    done
    [[ $c -eq 0 ]] && out+=( "${!i}" )
  done
  echo "${out[@]}"
}
export -f remove_duplicate

is_num() { #checks if all arguments are numbers
  local i res
  [[ -z $1 ]] && die "${FUNCNAME[0]}: Missing argument"
  for i in "$@"; do
    res=$(awk -v x="$i" 'BEGIN{ print ( x+0==x ) }') || die "${FUNCNAME[0]}: awk failed"
    [[ $res -eq 1 ]] || return 1
    unset res
  done
  return 0
}
export -f is_num

get_stepname() { #get the dir name of a certain step number (1st argument)
  local name
  [[ -n $1 ]] || die "${FUNCNAME[0]}: Missing argument"
  if [[ $1 = "--trunc" ]]; then
    echo "step_"
    return 0
  fi
  is_int "${1}" || die "${FUNCNAME[0]}: needs a int as argument, but got $1"
  name="$(printf step_%03i "$1")"
  [[ -z $name ]] && die "${FUNCNAME[0]}: Could not get stepname"
  echo "$name"
}
export -f get_stepname

update_stepnames(){ #updated the current working step to a certain number (1st argument)
  local thisstep laststep nr
  [[ -n $1 ]] || die "${FUNCNAME[0]}: Missing argument"
  nr="$1"
  is_int "$nr" || die "${FUNCNAME[0]}: needs a int as argument, but got $nr"
  [[ -z $CSG_MAINDIR ]] && die "${FUNCNAME[0]}: CSG_MAINDIR is undefined"
  [[ -d $CSG_MAINDIR ]] || die "${FUNCNAME[0]}: $CSG_MAINDIR is not dir"
  thisstep="$(get_stepname $nr)"
  export CSG_THISSTEP="$CSG_MAINDIR/$thisstep"
  if [[ $nr -gt 0 ]]; then
    laststep="$(get_stepname $((nr-1)) )"
    export CSG_LASTSTEP="$CSG_MAINDIR/$laststep"
  fi
}
export -f update_stepnames

get_current_step_dir() { #print the directory of the current step
  [[ -z $CSG_THISSTEP ]] && die "${FUNCNAME[0]}: \$CSG_THISSTEP is undefined (when calling from csg_call export it yourself)"
  if [[ $1 = "--no-check" ]]; then
    :
  else
    [[ -d $CSG_THISSTEP ]] || die "${FUNCNAME[0]}: $CSG_THISSTEP is not dir"
  fi
  echo "$CSG_THISSTEP"

}
export -f get_current_step_dir

get_last_step_dir() { #print the directory of the last step
  [[ -z $CSG_LASTSTEP ]] && die "${FUNCNAME[0]}: CSG_LASTSTEP is undefined  (when calling from csg_call export it yourself)"
  [[ -d $CSG_LASTSTEP ]] || die "${FUNCNAME[0]}: $CSG_LASTSTEP is not dir"
  echo "$CSG_LASTSTEP"
}
export -f get_last_step_dir

get_main_dir() { #print the main directory
  [[ -z $CSG_MAINDIR ]] && die "${FUNCNAME[0]}: CSG_MAINDIR is defined"
  [[ -d $CSG_MAINDIR ]] || die "${FUNCNAME[0]}: $CSG_MAINDIR is not dir"
  echo "$CSG_MAINDIR"
}
export -f get_main_dir

get_current_step_nr() { #print the main directory
  local name nr
  name=$(get_current_step_dir)
  nr=$(get_step_nr $name)
  echo "$nr"
}
export -f get_current_step_nr

get_step_nr() { #print the number of a certain step directory (1st argument)
  local nr trunc
  trunc=$(get_stepname --trunc)
  [[ -n $1 ]] || die "${FUNCNAME[0]}: Missing argument"
  nr=${1##*/}
  nr=${nr#$trunc}
  #convert to base 10 and cut leading zeros
  nr=$((10#$nr))
  is_int "$nr" || die "${FUNCNAME[0]}: Could not fetch step nr, got $nr"
  echo "$nr"
}
export -f get_step_nr

cp_from_main_dir() { #copy something from the main directory
  critical pushd "$(get_main_dir)"
  if [[ $1 = "--rename" ]]; then
    shift
    [[ $# -eq 2 && -n $1 && -n $2 ]] || die "${FUNCNAME[0]}: with --rename option has to be called with exactly 2 (non-empty) arguments"
    echo "cp_from_main_dir: '$1' to '$2'"
    critical cp $1 "$(dirs -l +1)/$2"
  else
    echo "cp_from_main_dir: '$@'"
    critical cp $@ "$(dirs -l +1)"
  fi
  critical popd
}
export -f cp_from_main_dir

cp_from_last_step() { #copy something from the last step
  if [[ $1 = "--rename" ]]; then
    shift
    [[ $# -eq 2 && -n $1 && -n $2 ]] || die "${FUNCNAME[0]}: with --rename option has to be called with exactly 2 (non-empty) arguments"
    echo "cp_from_last_step: '$1' to '$2'"
    critical pushd "$(get_last_step_dir)"
    critical cp $1 "$(dirs -l +1)/$2"
    critical popd
  else
    echo "cp_from_last_step: '$@'"
    critical pushd "$(get_last_step_dir)"
    critical cp $@ "$(dirs -l +1)"
    critical popd
  fi
}
export -f cp_from_last_step

get_time() { #gives back current time in sec from 1970
  date +%s || die "${FUNCNAME[0]}:  date +%s failed"
}
export -f get_time

get_number_tasks() { #get the number of possible tasks from the xml file or determine it automatically under some systems
  local tasks
  tasks="$(csg_get_property cg.inverse.simulation.tasks)"
  [[ $tasks = "auto" ]] && tasks=0
  is_int "$tasks" || die "${FUNCNAME[0]}: cg.inverse.simulation.tasks needs to be a number or 'auto', but I got $(csg_get_property cg.inverse.simulation.tasks)"
  if [[ $tasks -eq 0 ]]; then #auto-detect
    if [[ -r /proc/cpuinfo ]]; then #linux
      tasks=$(sed -n '/processor/p' /proc/cpuinfo | sed -n '$=')
    elif [[ -x /usr/sbin/sysctl ]]; then #mac os
      tasks=$(/usr/sbin/sysctl -n hw.ncpu)
    elif [[ -x /usr/sbin/lsdev ]]; then #AIX
      tasks=$(/usr/sbin/lsdev | sed -n '/Processor/p' | sed -n '$=')
    fi
    is_int "${tasks}" || tasks=1 #failback in case we got non-int
  fi
  if [[ ${CSG_NUM_THREADS} ]]; then
    is_int "${CSG_NUM_THREADS}" || die "${FUNCNAME[0]}: value of CSG_NUM_THREADS needs to be a number, but I got ${CSG_NUM_THREADS}"
    msg --color blue --to-stderr "${FUNCNAME[0]}: Overwriting cg.inverse.simulation.tasks with '${CSG_NUM_THREADS}'"
    tasks="${CSG_NUM_THREADS}"
  fi
  echo "$tasks"
}
export -f get_number_tasks

get_table_comment() { #get comment lines from a table and add common information, which include the git id and other information
  local version co
  [[ -n "$(type -p csg_call)" ]] || die "${FUNCNAME[0]}: Could not find csg_call"
  version="$(csg_call --version)" || die "${FUNCNAME[0]}: csg_call --version failed"
  echo "Created on $(date) by $USER@$HOSTNAME"
  echo "called from $version" | sed "s/csg_call/${0##*/}/"
  [[ -n "${CSGXMLFILE}" ]] && echo "settings file: '$(globalize_file "${CSGXMLFILE}")'"
  echo "working directory: $PWD"
  if [[ -f $1 ]]; then 
    co=$(sed -n 's/^[#@][[:space:]]*//p' "$1") || die "${FUNCNAME[0]}: sed failed"
    [[ -n $co ]] && echo "Comments from $(globalize_file $1):\n$co"
  fi
}
export -f get_table_comment

csg_inverse_clean() { #clean out the main directory 
  local i files log t
  [[ -n $1 ]] && t="$1" || t="30"
  log="$(csg_get_property cg.inverse.log_file 2>/dev/null)"
  echo -e "So, you want to clean?\n"
  echo "I will remove:"
  files="$(ls -d done ${log} $(get_stepname --trunc)* *~ 2>/dev/null)"
  if [[ -z $files ]]; then
    echo "Nothing to clean"
  else
    msg --color red $files
    msg --color blue "\nCTRL-C to stop it"
    for ((i=$t;i>0;i--)); do
      echo -n "$i "
      sleep 1
    done
    rm -rf $files
    msg --color green "\n\nDone, hope you are happy now"
  fi
}
export -f csg_inverse_clean

check_path_variable() { #check if a variable contains only valid paths
  local old_IFS dir
  [[ -z $1 ]] && die "${FUNCNAME[0]}: Missing argument"
  for var in "$@"; do
    [[ -z $var ]] && continue
    old_IFS="$IFS"
    IFS=":"
    for dir in ${!var}; do
      [[ -z $dir ]] && continue
      [[ $dir = *votca* ]] || continue #to many error otherwise
      [[ -d $dir ]] || die "${FUNCNAME[0]}: $dir from variable $var is not a directory"
    done
    IFS="$old_IFS"
  done
}
export -f check_path_variable

add_to_csgshare() { #added an directory to the csg internal search directories
  local dir end="no"
  [[ $1 = "--at-the-end" ]] && end="yes" && shift
  [[ -z $1 ]] && die "${FUNCNAME[0]}: Missing argument"
  for dirlist in "$@"; do
    old_IFS="$IFS"
    IFS=":"
    for dir in $dirlist; do
      #dir maybe contains $PWD or something
      eval dir="$dir"
      [[ -d $dir ]] || die "${FUNCNAME[0]}: Could not find scriptdir $dir"
      dir="$(globalize_dir "$dir")"
      if [[ $end = "yes" ]]; then
        export CSGSHARE="${CSGSHARE}${CSGSHARE:+:}$dir"
        export PERL5LIB="${PERL5LIB}${PERL5LIB:+:}$dir"
        export PYTHONPATH="${PYTHONPATH}${PYTHONPATH:+:}$dir"
      else
        export CSGSHARE="$dir${CSGSHARE:+:}$CSGSHARE"
        export PERL5LIB="$dir${PERL5LIB:+:}$PERL5LIB"
        export PYTHONPATH="$dir${PYTHONPATH:+:}$PYTHONPATH"
      fi
    done
    IFS="$old_IFS"
  done
  check_path_variable CSGSHARE PERL5LIB PYTHONPATH
}
export -f add_to_csgshare

globalize_dir() { #convert a local directory to a global one
  [[ -z $1 ]] && die "${FUNCNAME[0]}: missing argument"
  [[ -d $1 ]] || die "${FUNCNAME[0]}: '$1' is not a dir"
  cd "$1"
  pwd
}
export -f globalize_dir

globalize_file() { #convert a local file name to a global one
  [[ -z $1 ]] && die "${FUNCNAME[0]}: missing argument"
  [[ -f $1 ]] || die "${FUNCNAME[0]}: '$1' is not a file"
  local dir
  [[ ${1%/*} = ${1} ]] && dir="." || dir="${1%/*}"
  echo "$(globalize_dir "$dir")/${1##*/}"
}
export -f globalize_file

source_function() { #source an extra function file
  local function_file
  [[ -n $1 ]] || die "${FUNCNAME[0]}: Missing argument"
  function_file=$(source_wrapper functions $1) || die "${FUNCNAME[0]}: source_wrapper functions $1 failed"
  source ${function_file} || die "${FUNCNAME[0]}: source ${function_file} failed"
}
export -f source_function

csg_banner() { #print a big banner
  local i l=0 list=()
  [[ -z $1 ]] && return 0
  for i in "$@"; do
    while [[ -n $i && -z ${i/*\\n*} ]]; do
      list[$l]="${i%%\\n*}"
      ((l++))
      i="${i#*\\n}"
    done
    list[$l]=$i
    ((l++))
  done

  l="1"
  for i in "${list[@]}"; do
    [[ ${#l} -lt ${#i} ]] && l="${i}"
  done

  echo "####${l//?/#}"
  echo "# ${l//?/ } #"
  for i in "${list[@]}"; do
    printf "# %-${#l}s #\n" "$i"
  done
  echo "# ${l//?/ } #"
  echo "####${l//?/#}"
}
export -f csg_banner

csg_calc() { #simple calculator, a + b, ...
  local res ret=0 err="1e-2"
  [[ -z $1 || -z $2 || -z $3 ]] && die "${FUNCNAME[0]}: Needs 3 arguments, but got '$*'"
  is_num "$1" || die "${FUNCNAME[0]}: First argument of csg_calc should be a number, but got '$1'"
  is_num "$3" || die "${FUNCNAME[0]}: Third argument of csg_calc should be a number, but got '$3'"
  [[ -n "$(type -p awk)" ]] || die "${FUNCNAME[0]}: Could not find awk"
  #we use awk -v because then " 1 " or "1\n" is equal to 1
  case "$2" in
    "+"|"-"|'*'|"/"|"^")
       res="$(awk -v x="$1" -v y="$3" "BEGIN{print ( x $2 y ) }")" || die "${FUNCNAME[0]}: awk -v x='$1' -v y='$3' 'BEGIN{print ( x $2 y ) }' failed"
       true;;
    '>'|'<' )
       res="$(awk -v x="$1" -v y="$3" "BEGIN{print ( x $2 y )}")" || die "${FUNCNAME[0]}: awk -v x='$1' -v y='$3' 'BEGIN{print ( x $2 y )}' failed"
       #awk return 1 for true and 0 for false, shell exit codes are the other way around
       ret="$((1-$res))"
       #return value matters
       res=""
       true;;
    "="|"==")
       #this is really tricky... case x=0,y=0 is catched by (x==y) after that |x-y|/max(|x|,|y|) will work expect for x,y beginng close to zero
       res="$(awk -v x="$1" -v y="$3" -v e="$err" \
       'function max(x,y){return (x>y)?x:y;} function abs(x){return (x<0)?-x:x;} BEGIN{if (x==y){print 1;}else{if (abs(x-y)<e){print 1;}else{ print ( abs(x-y)/max(abs(x),abs(y)) < e );}}}')" \
	 || die "${FUNCNAME[0]}: awk for =/== failed"
       #awk return 1 for true and 0 for false, shell exit codes are the other way around
       ret="$((1-$res))"
       #return value matters
       res=""
       true;;
    *)
       die "${FUNCNAME[0]}: unknow operation" 
       true;;
  esac
  [[ -n $res ]] && echo "$res"
  return $ret
}
export -f csg_calc

show_csg_tables() { #show all concatinated csg tables
  local old_IFS dir
  old_IFS="$IFS"
  IFS=":"
  echo "#The order in which scripts get called"
  echo "#CSGSHARE is $CSGSHARE"
  for dir in ${CSGSHARE}; do
    [[ -f $dir/csg_table ]] || continue
    echo "#From: $dir/csg_table"
    #remove comments and empty line, trim begin and end, tab to spaces
    sed -e '/^#/d' -e '/^[[:space:]]*$/d' -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e 's/[[:space:]]\+/ /g' "$dir/csg_table"
  done
  IFS="$old_IFS"
}
export -f show_csg_tables

get_command_from_csg_tables() { #print the name of script belonging to certain tags (1st, 2nd argument)
  [[ -z $1 || -z $2 ]] && die "${FUNCNAME[0]}: Needs two tags"
  show_csg_tables | \
    sed -e '/^#/d' | \
    sed -n "s/^$1 $2 \(.*\)$/\1/p" | \
    sed -n '1p'
}
export -f get_command_from_csg_tables

source_wrapper() { #print the full name of a script belonging to two tags (1st, 2nd argument)
  [[ -z $1 || -z $2 ]] && die "${FUNCNAME[0]}: Needs two tags"
  local cmd script
  if [[ $1 = "function" ]]; then
    [[ $(type -t "$2") = "function" ]] || die "${FUNCNAME[0]}: could not find any function called '$2' (when calling from csg_call you might need to add --simprog option or set cg.inverse.program in the xml file)"
    echo "$2"
  else
    cmd=$(get_command_from_csg_tables "$1" "$2") || die "${FUNCNAME[0]}: get_command_from_csg_tables '$1' '$2' failed"
    [[ -z $cmd ]] && die "${FUNCNAME[0]}: Could not get any script from tags '$1' '$2'"
    #cmd might contain option after the script name
    script="${cmd%% *}"
    real_script="$(find_in_csgshare "$script")"
    echo "${cmd/${script}/${real_script}}"
  fi
}
export -f source_wrapper

find_in_csgshare() { #find a script in csg script search path
  [[ -z $1 ]] && die "${FUNCNAME[0]}: missing argument"
  #global path
  if [[ -z ${1##/*} ]]; then
    [[ -f $1 ]] || die "${FUNCNAME[0]}: $1 is a script with global path, but was not found"
    echo "$1" && return
  fi
  local old_IFS dir
  old_IFS="$IFS"
  IFS=":"
  for dir in ${CSGSHARE}; do
    [[ -f $dir/$1 ]] && break
  done
  IFS="$old_IFS"
  [[ -f $dir/$1 ]] && echo "$dir/$1" && return
  die "${FUNCNAME[0]}: Could not find script $1 in $CSGSHARE"
}
export -f find_in_csgshare

if [ -z "$(type -p mktemp)" ]; then
  #do not document this
  mktemp() {
    [[ $1 = "-u" ]] && shift
    [[ -z $1 ]] && die "${FUNCNAME[0]}: missing argument"
    [[ -z ${1##*X} ]] || die "${FUNCNAME[0]}: argument has to end at least with X"
    local end trunc i l tmp newend
    end=${1##*[^X]}
    trunc=${1%${end}}
    l=${end//[^X]}
    l=${#l}
    while true; do
      newend="$end"
      for ((i=0;i<$l;i++)); do
        newend="${newend/X/${RANDOM:0:1}}"
      done
      tmp="${trunc}${newend}"
      [[ -f $tmp ]] || break
    done
    echo "$tmp"
  }
  export -f mktemp
fi

enable_logging() { #enables the logging to a certain file (1st argument) or the logfile taken from the xml file
  local log
  if [[ -z $1 ]]; then
    log="$(csg_get_property cg.inverse.log_file 2>/dev/null)"
  else
    log="$1"
  fi
  log="${PWD}/${log##*/}"
  export CSGLOG="$log"
  if [[ -f $CSGLOG ]]; then
    exec 3>&1 4>&2 >> "$CSGLOG" 2>&1
    echo -e "\n\n#################################"
    echo "# Appending to existing logfile #"
    echo -e "#################################\n\n"
    msg --color blue "Appending to existing logfile ${CSGLOG##*/}"
  else
    exec 3>&1 4>&2 >> "$CSGLOG" 2>&1
    msg "For a more verbose log see: ${CSGLOG##*/}"
  fi
}
export -f enable_logging

get_restart_file() { #print the name of the restart file to use
  local file
  file="$(csg_get_property cg.inverse.restart_file)"
  [[ -z ${file/*\/*} ]] && die "${FUNCNAME[0]}: cg.inverse.restart_file has to be a local file with slash '/'"
  echo "$file"
}
export -f get_restart_file

check_for_obsolete_xml_options() { #check xml file for obsolete options
  local i
  for i in cg.inverse.mpi.tasks cg.inverse.mpi.cmd cg.inverse.parallel.tasks cg.inverse.parallel.cmd \
    cg.inverse.gromacs.mdrun.bin cg.inverse.espresso.bin cg.inverse.scriptdir cg.inverse.gromacs.grompp.topol \
    cg.inverse.gromacs.grompp.index cg.inverse.gromacs.g_rdf.topol cg.inverse.convergence_check \
    cg.inverse.convergence_check_options.name_glob cg.inverse.convergence_check_options.limit \
    cg.inverse.espresso.table_end cg.inverse.gromacs.traj_type cg.inverse.gromacs.topol_out \
    cg.inverse.espresso.blockfile cg.inverse.espresso.blockfile_out cg.inverse.espresso.n_steps \
    cg.inverse.espresso.exclusions cg.inverse.espresso.debug cg.inverse.espresso.n_snapshots \
    cg.non-bonded.inverse.espresso.index1 cg.non-bonded.inverse.espresso.index2 cg.inverse.espresso.success \
    cg.inverse.espresso.scriptdir cg.non-bonded.inverse.post_update_options.kbibi.type \
    cg.inverse.imc.numpy.bin cg.inverse.imc.octave.bin cg.inverse.imc.matlab.bin cg.inverse.imc.solver \
    cg.non-bonded.inverse.imc.reg \
    ; do
    [[ -z "$(csg_get_property --allow-empty $i)" ]] && continue #filter me away
    new=""
    case $i in
      cg.inverse.mpi.tasks|cg.inverse.parallel.tasks)
        new="cg.inverse.simulation.tasks";;
      cg.inverse.gromacs.mdrun.bin|cg.inverse.espresso.bin)
        new="${i/bin/command}";;
      cg.inverse.scriptdir)
        new="${i/dir/path}";;
      cg.inverse.gromacs.grompp.index)
        new="${i/.grompp}";;
      cg.inverse.gromacs.grompp.topol)
        new="cg.inverse.gromacs.topol_in";;
      cg.inverse.gromacs.g_rdf.topol)
        new="${i/g_}";;
      cg.inverse.gromacs.topol_out)
        new="${i/_out}";;
      cg.inverse.gromacs.traj_type)
        new="";;
      cg.inverse.convergence_check)
	new="${i}.type";;
      cg.inverse.convergence_check_options.limit)
        new="cg.inverse.convergence_check.limit";;
      cg.non-bonded.inverse.imc.reg)
        new="cg.inverse.imc.<group>.reg";;
    esac
    [[ -n $new ]] && new="has been renamed to $new" || new="has been removed"
    die "${FUNCNAME[0]}: The xml option $i $new\nPlease remove the obsolete options from the xmlfile"
  done
}
export -f check_for_obsolete_xml_options

command_not_found_handle() { #print and error message if a command or a function was not found
  die "Command/function $1 not found (when calling from csg_call you might need to add --simprog option or set cg.inverse.program in the xml file)"
}
export -f command_not_found_handle

#in bash4 this is not needed, but for older bash we add add a failback from most important simulation functions
for i in simulation_finish checkpoint_exist get_simulation_setting; do
  eval $i\(\) { command_not_found_handle $i\; }
  eval export -f $i 
done
unset i

