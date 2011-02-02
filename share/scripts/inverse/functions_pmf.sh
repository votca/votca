#! /bin/bash

run() {
  local logfile
  if [ "$1" = "--log" ]; then
    logfile="$2"
    shift 2
  else
    logfile="log_$1"
  fi
  echo Doing "$1 (Details see $logfile)"
  $* &> "$logfile" || die "$1 failed"
}

die() {
  echo "$*" >&2
  exit 1
}

get_pullgroup_nr() {
  local name appearance number confin
  [ -z "$2" ] && die "get_pullgroup_nr: Missing argument"
  name="$1"
  confin="$2"
  if [ -z "$3" ]; then
    appearance="1"
  else
    appearance="$3"
  fi
  number="$(grep "$name" ${confin} | sed -n "${appearance}p" | perl -pe "s/^.*$name\s*([0-9]+)\s+.*$/\$1/" )" || die "get_pullgroup_nr: failed"
  [ -z "$number" ] && die "get_pullgroup_nr: Could not get number"
  echo $number
}

get_from_mdp() {
  local res
  [[ -n "$2" ]] || die "get_from_mdp: Missing argument (what file)"
  [[ -f "$2" ]] || die "get_from_mdp: Could not read file '$2'"
  res=$(sed -n -e "s#[[:space:]]*$1[[:space:]]*=[[:space:]]*\(.*\)\$#\1#p" $2 | sed -e 's#;.*##' -e 's#[[:space:]]##g' ) || \
    die "get_from_mdp: sed failed"
  [[ -n "$res" ]] || die "get_from_mdp: could not fetch $1 from $2"
  echo "$res"
}

int_check(){
  [[ -n "$1" ]] || die "int_check: Missig argument"
  [[ -z "${1//[0-9]}" ]] && return 0
  return 1      
}

get_step_nr(){ 
  local number
  number=${PWD##*/}
  number=${number%%_*}
  int_check "$number" || die "get_step_nr: Could not fetch nr (number was '$number')"
  echo $number
}

get_step_name(){
  local nicenumber name
  [ -z "$1" ] && die "get_step_name: Missing argument"
  int_check "$1" || die "get_step_name: Argument should be a number"
  nicenumber=$(printf %02i $1)
  cd ..
  name=$(ls -d ${nicenumber}_*) || die "get_step_name: Could not find step nr $1"
  echo $name
}

[ "$fake_qstart" = "yes" ] && \
qstart() {
  local ffile
  while [ "$1" != "--" ] && [ -n "$1" ]; do
    [ "$1" = "-r" ] || ffile="$2"
    shift 
  done
  shift
  [ -n "$1" ] || die "fake_qstart: Missing argument"
  echo -e "qstart: running $* (details see log_qstart)"
  $* >> log_qstart 2>&1 || die "qstart: '$*' failed"
  [ -n "$ffile" ] && touch $ffile
}

cp_from_last_step_to() {
  [ -d "$last_step" ] || die "Last step dir $last_step not found"
  [ -z "$2" ] && die "cp_from_last_step_to: missing argument"
  cp $last_step/$1 $2
}

cp_from_last_step_here() {
  local i
  [ -z "$1" ] && die "cp_form_last_step: missing argument"
  for i in $@; do
    cp_from_last_step_to $i . 
  done
}
