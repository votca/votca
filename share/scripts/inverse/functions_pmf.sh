#! /bin/bash

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
