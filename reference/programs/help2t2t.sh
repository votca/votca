#! /bin/bash

die() {
  echo -e "$*" >&2
  exit 1
}

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing argument"
prog="${1}"
[[ -f $2 ]] || die "${0##*/}: Could not open '$2'"
helpmsg="$(<$2)"

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
#-examplelines (^*) -> - ``command`` descrition
#  descrition is optional and start after >=2 spaces
#-extra empty line before new section to close itemize
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
   -e 's/^\(Threading\|Allowed\|Trajectory\|Mapping\|Specific\) options:/\n&/'

