#! /bin/bash

die() {
  echo -e "$*" >&2
  exit 1
}

[ -z "$1" ] && die "Give the path to votca src as argument"
[ -d "$1" ] || die "Arg '$1' is not a dir"

find "$1" -name "*.cc" -o -name "*.h" \
          -exec grep boost {} \; | \
grep include | \
sed 's/^.*<\([^>]*\)>.*$/\1/'
