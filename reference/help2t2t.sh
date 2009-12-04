#! /bin/bash

die() {
  echo -e "$*" >&2
  exit 1
}

[ -z "$1" ] && die "${0##*/}: Missing argument"
prog="$1"
shift

[ -z "${CSGSHARE}" ] && die "${0##*/}: CSGSHARE not defined"
 
[ -z "$(type $prog)" ] && die "${0##*/}: $prog not found"

helpmsg="$($prog "$@" --help)" || die "${0##*/}: $prog $@ --help failed"

echo "$prog $@"
echo ${0##*/}
date
echo
echo '%!includeconf: config.t2t'
echo 

echo -e "++$prog++"
echo -e "$helpmsg" | sed \
   -e '1,3d' \
   -e 's/^[[:space:]]*//' \
   -e '/^-/s/ \{2\}/``  /' \
   -e '/^-/s/^/- ``/' \
   -e 's/^\(Usage:[[:space:]]*\)\(.*\)$/\1``\2``/' \
   -e 's/^\(Examples\|USES\|NEEDS|Usage\):/\n&/' \
   -e 's/^\* \(.*\)$/- ``\1``/' \

