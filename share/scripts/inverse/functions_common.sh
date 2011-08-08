#! /bin/bash
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

msg() { #echos a msg on the screen and send it to the logfile if logging is enabled
  local color colors=" blue cyan cyann green red purp "
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
    [[ -z $2 ]] && die "msg: missing argument after --color"
    [[ -n ${colors//* $2 *} ]] && die "msg: Unknown color ($colors allowed)"
    color="${!2}"
    shift 2
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

unset -f die
die () { #make the iterative frame work stopp
  local pid pids c
  msg --color red --to-stderr "$(csg_banner "ERROR:" "$@")"
  [[ -z $CSGLOG ]] || msg --color blue "For details see $CSGLOG"
  if [[ -n ${CSG_MASTER_PID} ]]; then
    #grabbing the pid group would be easier, but it would not work on AIX
    pid="$$"
    pids="$$"
    c=0
    #find the parent of pid until we reach CSG_MASTER_PID
    until [[ ${CSG_MASTER_PID} -eq $pid ]]; do
      #get the parent pid using BSD style due to AIX
      pid=$(ps -o ppid= -p "$pid" 2>/dev/null)
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
      echo "die: (called from $$)  CSG_MASTER_PID is $CSG_MASTER_PID" >&2
      echo "die: pids to kill: $pids" >&2
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
  script="$(source_wrapper $1 $2)" || die "cat_external: source_wrapper $1 $2 failed"
  cat "${script/ *}"
}
export -f cat_external

do_external() { #takes two tags, find the according script and excute it
  local script tags quiet="no" 
  [[ $1 = "-q" ]] && quiet="yes" && shift
  script="$(source_wrapper $1 $2)" || die "do_external: source_wrapper $1 $2 failed"
  tags="$1 $2"
  shift 2

  [[ $quiet = "no" ]] && echo "Running subscript '${script##*/} $*' (from tags $tags) dir ${script%/*}"
  if [[ -n $CSGDEBUG && -n "$(sed -n '1s@bash@XXX@p' "$script")" ]]; then
    bash -x $script "$@"
  elif [[ -n $CSGDEBUG && -n "$(sed -n '1s@perl@XXX@p' "$script")" ]]; then
    local perl_debug="$(mktemp perl_debug.XXX)" ret
    PERLDB_OPTS="NonStop=1 AutoTrace=1 frame=2 LineInfo=$perl_debug" perl -dS $script "$@"
    ret=$?
    cat "$perl_debug" 2>&1
    [[ $ret -eq 0 ]]
  else
    $script "$@"
  fi || die "do_external: subscript $script $* (from tags $tags) failed"
}
export -f do_external

critical() { #executes arguments as command and calls die if not succesful
  local quiet="no"
  [[ $1 = "-q" ]] && quiet="yes" && shift
  [[ -z $1 ]] && die "critical: missing argument"
  #print this message to stderr because $(critical something) is used very often
  [[ $quiet = "no" ]] && echo "Running critical command '$*'" >&2
   "$@" || die "critical: '$*' failed"
}
export -f critical

for_all (){ #do something for all interactions (1st argument)
  local bondtype name interactions quiet="no"
  [[ $1 = "-q" ]] && quiet="yes" && shift
  [[ -z $1 || -z $2 ]] && "for_all need at least two arguments"
  bondtype="$1"
  shift
  #check that type is bonded or non-bonded
  if [[ $bondtype != "non-bonded" ]]; then
    die  "for_all: Argmuent 1 '$bondtype' is not non-bonded"
  fi
  [[ $quiet = "no" ]] && echo "For all $bondtype" >&2
  check_for_duplicated_interactions
  interactions="$(csg_get_property cg.non-bonded.name)"
  for name in $interactions; do
    #print this message to stderr to avoid problem with $(for_all something)
    [[ $quiet = no ]] && echo "for_all: run '$*'" >&2
    #we need to use bash -c here to allow things like $(csg_get_interaction_property name) in arguments
    #write variable defines in the front is better, that export
    #no need to run unset afterwards
    bondtype="$bondtype" \
    bondname="$name" \
    bash -c "$*" || die "for_all: bash -c '$*' failed for bondname '$name'"
  done
}
export -f for_all

check_for_duplicated_interactions() { #checks for duplicated interactions
  local i j names=( $(csg_get_property cg.non-bonded.name) ) 
  for ((i=0;i<${#names[@]};i++)); do
    for ((j=i+1;j<${#names[@]};j++)); do
      [[ ${names[$i]} = ${names[$j]} ]] && die "for_all: the interaction name '${names[$i]}' appeared twice, this is not allowed"
    done
  done
}
export -f check_for_duplicated_interactions

csg_get_interaction_property () { #gets an interaction property from the xml file, should only be called from inside a for_all loop
  local ret allow_empty cmd
  if [[ $1 = "--allow-empty" ]]; then
    shift
    allow_empty="yes"
  else
    allow_empty="no"
  fi
  [[ -n $1 ]] || die "csg_get_interaction_property: Missing argument"
  [[ -n $bondtype ]] || die "csg_get_interaction_property: bondtype is undefined (when calling from csg_call set it by --ia-type option)"
  
  #bondtype is special -> dirty hack - removed whenever issue 13 is fixed
  #make these this case work even without name (called by csg_call)
  [[ $1 = "bondtype" ]] && echo "$bondtype" && return 0

  [[ -n $bondname ]] || die "csg_get_interaction_property: bondname is undefined (when calling from csg_call set it by --ia-name option)"
  #make these this case work even without xml file (called by csg_call)
  [[ $1 = "name" ]] && echo "$bondname" && return 0

  [[ -n $CSGXMLFILE ]] || die "csg_get_interaction_property: CSGXMLFILE is undefined (when calling from csg_call set it by --options option)"
  [[ -n "$(type -p csg_property)" ]] || die "csg_get_interaction_property: Could not find csg_property"
  cmd="csg_property --file $CSGXMLFILE --short --path cg.${bondtype} --filter name=$bondname --print $1"
  #the --filter option will make csg_property fail if $1 does not exist, don't stop if we have an default
  if ! ret="$($cmd)"; then
    [[ $allow_empty = "no" && -z $2 ]] && \
      die "csg_get_interaction_property:\n'$cmd'\nfailed geting '$1' for interaction '$bondname' with error msg:\n $ret\n and no default for $1"
    #ret has error message
    ret=""
  fi
  ret="${ret%%[[:space:]]}"
  ret="${ret##[[:space:]]}"
  [[ $allow_empty = no && -z $ret && -n $2 ]] && ret="$2"
  [[ $allow_empty = no && -z $ret ]] && die "csg_get_interaction_property: Could not get '$1' for interaction '$bondname'\nResult of '$cmd' was empty"
  echo "${ret}"
}
export -f csg_get_interaction_property

csg_get_property () { #get an property from the xml file
  local ret allow_empty cmd
  if [[ $1 = "--allow-empty" ]]; then
    shift
    allow_empty="yes"
  else
    allow_empty="no"
  fi
  [[ -n $1 ]] || die "csg_get_property: Missing argument"
  [[ -n $CSGXMLFILE ]] || die "csg_get_property: CSGXMLFILE is undefined (when calling from csg_call set it by --options option)"
  [[ -n "$(type -p csg_property)" ]] || die "csg_get_property: Could not find csg_property"
  cmd="csg_property --file $CSGXMLFILE --path ${1} --short --print ."
  #csg_property only fails if xml file is bad otherwise result is empty
  ret="$(critical -q $cmd)"
  ret="${ret%%[[:space:]]}"
  ret="${ret##[[:space:]]}"
  [[ -z $ret && -n $2 ]] && ret="$2"
  [[ $allow_empty = "no" && -z $ret ]] && die "csg_get_property: Could not get '$1'\nResult of '$cmd' was empty"
  echo "${ret}"
}
export -f csg_get_property

mark_done () { #mark a task (1st argument) as done in the restart file
  local file
  [[ -n $1 ]] || die "mark_done: Missing argument"
  file="$(get_restart_file)"
  is_done "$1" || echo "$1 done" >> "${file}"
}
export -f mark_done

is_done () { #checks if something is already do in the restart file
  local file
  [[ -n $1 ]] || die "is_done: Missing argument"
  file="$(get_restart_file)"
  [[ -f ${file} ]] || return 1
  [[ -n "$(sed -n "/^$1 done\$/p" ${file})" ]] && return 0
  return 1
}
export -f is_done

is_int() { #checks if all arguments are integers
  local i
  [[ -z $1 ]] && die "is_int: Missing argument"
  for i in "$@"; do
    [[ -n $i && -z ${i//[0-9]} ]] || return 1
  done
  return 0
}
export -f is_int

is_num() { #checks if all arguments are numbers
  local i res
  [[ -z $1 ]] && die "is_num: Missing argument"
  for i in "$@"; do
    res=$(awk -v x="$i" 'BEGIN{ print x+0==x; }')
    [[ $res -eq 1 ]] || return 1
    unset res
  done
  return 0
}
export -f is_num

get_stepname() { #get the dir name of a certain step number (1st argument)
  local name
  [[ -n $1 ]] || die "get_stepname: Missing argument"
  if [[ $1 = "--trunc" ]]; then
    echo "step_"
    return 0
  fi
  is_int "${1}" || die "get_stepname: needs a int as argument, but got $1"
  name="$(printf step_%03i "$1")"
  [[ -z $name ]] && die "get_stepname: Could not get stepname"
  echo "$name"
}
export -f get_stepname

update_stepnames(){ #updated the current working step to a certain number (1st argument)
  local thisstep laststep nr
  [[ -n $1 ]] || die "update_stepnames: Missing argument"
  nr="$1"
  is_int "$nr" || die "update_stepnames: needs a int as argument, but got $nr"
  [[ -z $CSG_MAINDIR ]] && die "update_stepnames: CSG_MAINDIR is undefined"
  [[ -d $CSG_MAINDIR ]] || die "update_stepnames: $CSG_MAINDIR is not dir"
  thisstep="$(get_stepname $nr)"
  export CSG_THISSTEP="$CSG_MAINDIR/$thisstep"
  if [[ $nr -gt 0 ]]; then
    laststep="$(get_stepname $((nr-1)) )"
    export CSG_LASTSTEP="$CSG_MAINDIR/$laststep"
  fi
}
export -f update_stepnames

get_current_step_dir() { #print the directory of the current step
  [[ -z $CSG_THISSTEP ]] && die "get_current_step_dir: \$CSG_THISSTEP is undefined (when calling from csg_call export it yourself)"
  if [[ $1 = "--no-check" ]]; then
    :
  else
    [[ -d $CSG_THISSTEP ]] || die "get_last_step_dir: $CSG_THISSTEP is not dir"
  fi
  echo "$CSG_THISSTEP"

}
export -f get_current_step_dir

get_last_step_dir() { #print the directory of the last step
  [[ -z $CSG_LASTSTEP ]] && die "get_last_step_dir: CSG_LASTSTEP is undefined  (when calling from csg_call export it yourself)"
  [[ -d $CSG_LASTSTEP ]] || die "get_last_step_dir: $CSG_LASTSTEP is not dir"
  echo "$CSG_LASTSTEP"
}
export -f get_last_step_dir

get_main_dir() { #print the main directory
  [[ -z $CSG_MAINDIR ]] && die "get_main_dir: CSG_MAINDIR is defined"
  [[ -d $CSG_MAINDIR ]] || die "update_stepnames: $CSG_MAINDIR is not dir"
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
  [[ -n $1 ]] || die "get_step_nr: Missing argument"
  nr=${1##*/}
  nr=${nr#$trunc}
  #convert to base 10 and cut leading zeros
  nr=$((10#$nr))
  is_int "$nr" || die "get_step_nr: Could not fetch step nr, got $nr"
  echo "$nr"
}
export -f get_step_nr

cp_from_main_dir() { #copy something from the main directory
  if [[ $1 = "--rename" ]]; then
    shift
    [[ $# -eq 2 && -n $1 && -n $2 ]] || die "cp_from_main_dir: with --rename option has to be called with exactly 2 (non-empty) arguments"
    echo "cp_from_main_dir: '$1' to '$2'"
    critical pushd "$(get_main_dir)"
    critical cp $1 "$(dirs -l +1)/$2"
    critical popd
  else
    echo "cp_from_main_dir: '$@'"
    critical pushd "$(get_main_dir)"
    critical cp $@ "$(dirs -l +1)"
    critical popd
  fi
}
export -f cp_from_main_dir

cp_from_last_step() { #copy something from the last step
  if [[ $1 = "--rename" ]]; then
    shift
    [[ $# -eq 2 && -n $1 && -n $2 ]] || die "cp_from_last_step: with --rename option has to be called with exactly 2 (non-empty) arguments"
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

get_time() {
  date +%s || die "get_time:  date +%s failed"
}
export -f get_time

get_number_tasks() { #get the number of possible tasks from the xml file or determine it automatically under linux
  local tasks
  tasks="$(csg_get_property cg.inverse.simulation.tasks "auto")"
  [[ $tasks = "auto" ]] && tasks=0
  is_int "$tasks" || die "get_number_tasks: cg.inverse.parallel.tasks needs to be a number or 'auto', but I got $tasks"
  #this only work for linux
  if [[ $tasks -eq 0 && -r /proc/cpuinfo ]]; then
    tasks=$(sed -n '/processor/p' /proc/cpuinfo | sed -n '$=')
    [[ -z ${tasks//[0-9]} ]] || tasks=1
  fi
  [[ $tasks -le 1 ]] && tasks=1
  echo "$tasks"
}
export -f get_number_tasks

get_table_comment() { #get comment lines from a table and add common information, which include the hgid and other information
  local version co
  [[ -n "$(type -p csg_call)" ]] || die "get_defaults_comment: Could not find csg_version"
  version="$(csg_call --version)" || die "get_defaults_comment: csg_call --version failed"
  echo "Created on $(date) by $USER@$HOSTNAME"
  echo "called from $version" | sed "s/csg_call/${0##*/}/"
  [[ -n ${CSGXMLFILE} ]] && echo "settings file: $(globalize_file $CSGXMLFILE)"
  echo "working directory: $PWD"
  if [[ -f $1 ]]; then 
    co=$(sed -n 's/^[#@][[:space:]]*//p' "$1") || die "get_table_comment: sed failed"
    [[ -n $co ]] && echo "Comments from $(globalize_file $1):\n$co"
  fi
}
export -f get_table_comment

csg_inverse_clean() { #clean out the main directory 
  local i files log t
  [[ -n $1 ]] && t="$1" || t="30"
  log="$(csg_get_property cg.inverse.log_file "inverse.log")"
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
  [[ -z $1 ]] && die "check_path_variable: Missing argument"
  for var in "$@"; do
    [[ -z $var ]] && continue
    old_IFS="$IFS"
    IFS=":"
    for dir in ${!var}; do
      [[ -z $dir ]] && continue
      [[ -d $dir ]] || die "check_path_variable: $dir from variable $var is not a directory"
    done
    IFS="$old_IFS"
  done
}
export -f check_path_variable

add_to_csgshare() { #added an directory to the csg internal search directories
  local dir end="no"
  [[ $1 = "--at-the-end" ]] && end="yes" && shift
  [[ -z $1 ]] && die "add_to_csgshare: Missing argument"
  for dirlist in "$@"; do
    old_IFS="$IFS"
    IFS=":"
    for dir in $dirlist; do
      #dir maybe contains $PWD or something
      eval dir="$dir"
      [[ -d $dir ]] || die "add_to_csgshare: Could not find scriptdir $dir"
      dir="$(globalize_dir "$dir")"
      if [[ $end = "yes" ]]; then
        export CSGSHARE="${CSGSHARE}${CSGSHARE:+:}$dir"
        export PERL5LIB="${PERL5LIB}${PERL5LIB:+:}$dir"
      else
        export CSGSHARE="$dir${CSGSHARE:+:}$CSGSHARE"
        export PERL5LIB="$dir${PERL5LIB:+:}$PERL5LIB"
      fi
    done
    IFS="$old_IFS"
  done
  check_path_variable CSGSHARE PERL5LIB
}
export -f add_to_csgshare

globalize_dir() { #convert a local directory to a global one
  [[ -z $1 ]] && die "globalize_dir: missing argument"
  [[ -d $1 ]] || die "globalize_dir: '$1' is not a dir"
  cd "$1"
  pwd
}
export -f globalize_dir

globalize_file() { #convert a local file name to a global one
  [[ -z $1 ]] && die "globalize_file: missing argument"
  [[ -f $1 ]] || die "globalize_file: '$1' is not a file"
  local dir
  [[ ${1%/*} = ${1} ]] && dir="." || dir="${1%/*}"
  echo "$(globalize_dir "$dir")/${1##*/}"
}
export -f globalize_file

source_function() { #source an extra function file
  local function_file
  [[ -n $1 ]] || die "source_function: Missing argument"
  function_file=$(source_wrapper functions $1) || die "source_function: source_wrapper functions $1 failed"
  source ${function_file} || die "source_function: source ${function_file} failed"
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
  [[ -z $1 || -z $2 || -z $3 ]] && die "csg_calc: Needs 3 arguments, but got '$*'"
  is_num "$1" || die "csg_calc: First argument of csg_calc should be a number, but got '$1'"
  is_num "$3" || die "csg_calc: Third argument of csg_calc should be a number, but got '$3'"
  [[ -n "$(type -p awk)" ]] || die "csg_calc: Could not find awk"
  #we use awk -v because then " 1 " or "1\n" is equal to 1
  case "$2" in
    "+"|"-"|'*'|"/"|"**")
       res="$(awk -v x="$1" -v y="$3" "BEGIN{print x $2 y}")" || die "csg_calc: awk -v x='$1' -v y='$3' 'BEGIN{print x $2 y}' failed"
       true;;
    '>'|'<' )
       res="$(awk -v x="$1" -v y="$3" "BEGIN{print ( x $2 y )}")" || die "csg_calc: awk -v x='$1' -v y='$3' 'BEGIN{print ( x $2 y )}' failed"
       #awk return 1 for true and 0 for false, shell exit codes are the other way around
       ret="$((1-$res))"
       #return value matters
       res=""
       true;;
    "="|"==")
       #we expect that x and y are close together
       res="$(awk -v x="$1" -v y="$3" "BEGIN{print ( (sqrt((x-y)/x)**2) < $err )}")" || die "csg_calc: awk -v x='$1' -v y='$3' 'BEGIN{print ( (sqrt((x-y)/x)**2) < $err )}' failed"
       #awk return 1 for true and 0 for false, shell exit codes are the other way around
       ret="$((1-$res))"
       #return value matters
       res=""
       true;;
    *)
       die "csg_calc: unknow operation" 
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
  [[ -z $1 || -z $2 ]] && die "get_command_from_csg_tables: Needs two tags"
  show_csg_tables | \
    sed -e '/^#/d' | \
    sed -n "s/^$1 $2 \(.*\)$/\1/p" | \
    sed -n '1p'
}
export -f get_command_from_csg_tables

source_wrapper() { #print the full name of a script belonging to two tags (1st, 2nd argument)
  [[ -z $1 || -z $2 ]] && die "source_wrapper: Needs two tags"
  local cmd script
  cmd=$(get_command_from_csg_tables "$1" "$2") || die
  [[ -z $cmd ]] && die "source_wrapper: Could not get any script from tags '$1' '$2'"
  script="${cmd/* }"
  real_script="$(find_in_csgshare "$script")"
  echo "${cmd/${script}/${real_script}}"
}
export -f source_wrapper

find_in_csgshare() { #find a script in csg script search path
  [[ -z $1 ]] && die "find_in_csgshare: missing argument"
  #global path
  if [[ -z ${1##/*} ]]; then
    [[ -f $1 ]] || die "find_in_csgshare: $1 is a script with global path, but was not found"
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
  die "find_in_csgshare: Could not find script $1 in $CSGSHARE"
}
export -f find_in_csgshare

if [ -z "$(type -p mktemp)" ]; then
  #do not document this
  mktemp() {
    [[ -z $1 ]] && die "mktemp: missing argument"
    [[ -z ${1##*X} ]] || die "mktemp: argument has to end at least with X"
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
    log="$(csg_get_property cg.inverse.log_file "inverse.log")"
  else
    log="$1"
  fi
  log="${PWD}/${log##*/}"
  export CSGLOG="$log"
  if [[ -f $CSGLOG ]]; then
    exec 3>&1 4>&2 >> "$CSGLOG" 2>&1
    echo "\n\n#################################"
    echo "# Appending to existing logfile #"
    echo "#################################\n\n"
    msg --color blue "Appending to existing logfile ${CSGLOG##*/}"
  else
    exec 3>&1 4>&2 >> "$CSGLOG" 2>&1
    msg "For a more verbose log see: ${CSGLOG##*/}"
  fi
}
export -f enable_logging

get_restart_file() { #print the name of the restart file to use
  local file
  file="$(csg_get_property cg.inverse.restart_file "restart_points.log")"
  [[ -z ${file/*\/*} ]] && die "get_restart_file: cg.inverse.restart_file has to be a local file with slash '/'"
  echo "$file"
}
export -f get_restart_file

check_for_obsolete_xml_options() { #check xml file for obsolete options
  local i
  for i in cg.inverse.mpi.tasks cg.inverse.mpi.cmd cg.inverse.parallel.tasks cg.inverse.parallel.cmd \
    cg.inverse.gromacs.mdrun.bin cg.inverse.espresso.bin cg.inverse.scriptdir; do
    [[ -z "$(csg_get_property --allow-empty $i)" ]] && continue #filter me away
    case $i in
      cg.inverse.parallel.cmd|cg.inverse.mpi.cmd)
        new="";;
      cg.inverse.mpi.tasks|cg.inverse.parallel.tasks)
        new="cg.inverse.simulation.tasks";;
      cg.inverse.gromacs.mdrun.bin|cg.inverse.espresso.bin)
        new="${i/bin/command}";;
      cg.inverse.scriptdir)
        new="${i/dir/path}";;
      *)
        die "check_for_obsolete_xml_options: Unknown new name for obsolete xml option '$i'";;
    esac
    [[ -n $new ]] && new="has been renamed to $new" || new="has been removed"
    die "The xml option $i $new\nPlease remove the obsolete options from the xmlfile"
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

