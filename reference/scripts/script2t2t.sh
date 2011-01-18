#! /bin/bash

die() {
  echo -e "$*" >&2
  exit 1
}

[ -z "$1" ] && die "${0##*/}: Missing argument"
script="$1"
shift

[ -z "$(type csg_call)" ] && die "${0##*/}: csg_call not found"
helpmsg="$(csg_call --direct "$script" --help)" || die "${0##*/}: csg_call --direct failed"

echo "$script"
echo ${0##*/}
date
echo
echo '%!includeconf: config.t2t'
echo 
echo -e "++$script++"
echo -e "label($script)"


#Here comes the magic:
#-header in trash
#-space at begining and end are removed
#-optionlines (-[^ ]) -> - ``--option`` text
#-usageline -> usage: ``code``
#-convert NEEDS and OPTINAL in link
#-examplelines (^*) -> - ``line``
#-add cg.interaction in links
#-extra empty line before new section to close itemize
echo -e "$helpmsg" | sed \
   -e '1,/^please submit/d' \
   -e "/^$script/d" \
   -e 's/^[[:space:]]*//' \
   -e 's/[[:space:]]*$//' \
   -e '/^-[^ ]/s/ \{2\}/`` /' \
   -e '/^-[^ ].*``/s/^/- ``/' \
   -e 's/^\(Usage:[[:space:]]*\)\(.*\)$/\1``\2``/' \
   -e '/^\(NEEDS\|OPTIONAL\):/s/\([[:space:]]\)\([^[:space:]]\+\)/\1link(\2)(\2)/g' \
   -e 's/link(\([^c][^g][^)]*\))(\([^)]*\))/link(cg.interaction.\1)(\2)/g' \
   -e '/^\* /s/\( \{2\}\|$\)/`` /' \
   -e '/^\*.*``/s/^\*[[:space:]]*/- ``/' \
   -e 's/^\(Examples\|USES\|NEEDS\|Usage\|PROVIDES\|OPTIONAL\):/\n&/' \
   -e 's/^\(Allowed\|Trajectory\|Specific\) options:/\n&/' \


exit
#do not parse xml options of functions files
[ -z "${script##function_*}" ] && exit 0

content="$(csg_call --cat "$script")" || die "${0##*/}: csg_call --cat failed"
helpmsg="$(echo -e "$content" | sed -n -e '/^USES/d' -e '/csg_get_property/p')" || die "${0##*/}: sed failed"
#shell sciprts
if [ -n "$helpmsg" ] && [ -z "${script%%*.sh}" ]; then
  echo "Used xml properties:"
  #trick to manage multiple csg_get_property per line
  echo -e "$helpmsg" | sed -e 's/([[:space:]]\(csg_get_property\)/\n(\1/g'

