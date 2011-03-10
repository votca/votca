#! /bin/bash

die() {
  echo -e "$*" >&2
  exit 1
}

[ -z "$1" ] && die "${0##*/}: Missing argument"
prog="$1"
shift

[ -z "$(type $prog)" ] && die "${0##*/}: $prog not found"

helpmsg="$($prog "$@" --help)" || die "${0##*/}: $prog --help failed"

echo "$prog"
echo ${0##*/}
date
echo
echo '%!includeconf: config.t2t'
echo 
echo -e "++$prog++"
echo -e "label($prog)"

#Here comes the magic:
#-header in trash
#-remove prog/tools/gromacs version line
#-space at begining and end are removed
#-optionlines (-[^ ]) -> - ``--option`` text
#-usageline -> usage: ``code``
#-extra empty line before new section to close itemize
#-examplelines (^*) -> - ``line``
#-add cg.interaction in links
echo -e "$helpmsg" | sed \
   -e '1,/^please submit/d' \
   -e "/^\($prog\|votca_tools\|gromacs\)/d" \
   -e 's/^[[:space:]]*//' \
   -e 's/[[:space:]]*$//' \
   -e '/^-[^ ]/s/ \{2\}/`` /' \
   -e '/^-[^ ].*``/s/^/- ``/' \
   -e 's/^\(Usage:[[:space:]]*\)\(.*\)$/\1``\2``/' \
   -e '/^\* /s/\( \{2\}\|$\)/`` /' \
   -e '/^\*.*``/s/^\*[[:space:]]*/- ``/' \
   -e 's/^\(Examples\|Usage\):/\n&/' \
   -e 's/^\(Allowed\|Trajectory\|Specific\) options:/\n&/' \

