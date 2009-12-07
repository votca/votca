#! /bin/bash

die() {
  echo -e "$*" >&2
  exit 1
}

[ -z "$1" ] && die "${0##*/}: Missing argument"
prog="$1"
shift

if [ "$1" = "--direct" ]; then
  [ -z "$2" ] && die "${0##*/}: Missing argument after --direct"
  progname="$2"
else
  progname="$prog"
fi

[ -z "${CSGSHARE}" ] && die "${0##*/}: CSGSHARE not defined"
 
[ -z "$(type $prog)" ] && die "${0##*/}: $prog not found"

helpmsg="$($prog "$@" --help)" || die "${0##*/}: $progname --help failed"

echo "$progname"
echo ${0##*/}
date
echo
echo '%!includeconf: config.t2t'
echo 
echo -e "++$progname++"
echo -e "label($progname)"

#Here comes the magic:
#-header is trash
#-space at begining and end are removed
#-optionlines (-[^ ]) -> - ``--option`` text
#-usageline -> usage: ``code``
#-extra empty line before new section to close itemize
#-examplelines (^*) -> - ``line``
#-
echo -e "$helpmsg" | sed \
   -e '1,3d' \
   -e 's/^[[:space:]]*//' \
   -e 's/[[:space:]]*$//' \
   -e '/^-[^ ]/s/ \{2\}/`` /' \
   -e '/^-[^ ].*``/s/^/- ``/' \
   -e 's/^\(Usage:[[:space:]]*\)\(.*\)$/\1``\2``/' \
   -e 's/^\(Examples\|USES\|NEEDS\|Usage\|PROVIDES\|OPTIONAL\):/\n&/' \
   -e 's/^\(Allowed\|Trajectory\|Specific\) options:/\n&/' \
   -e '/^\* /s/\( \{2\}\|$\)/`` /' \
   -e '/^\*.*``/s/^\*[[:space:]]*/- ``/' \

