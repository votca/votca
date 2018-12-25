#! /bin/bash

die() {
  echo -e "$*" >&2
  exit 1
}

[[ -z $1 ]] && die "${0##*/}: Missing argument"
[[ -f $1 ]] || die "${0##*/}: Could not open '$1'"

echo "csg_table"
echo "make"
date
echo
echo '%!includeconf: config.t2t'
echo
echo "|| Key1 | Key2 | Scriptname"
cat "$1" | \
sed -e '/^#/d' | \
awk '{printf "| %s | %s | ref(%s)(%s) |\n",$1,$2,$3,$3}'
