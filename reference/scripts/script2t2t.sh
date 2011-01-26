#! /bin/bash

die() {
  echo -e "$*" >&2
  exit 1
}

[ -z "$1" ] && die "${0##*/}: Missing argument"
script="$1"
shift

[ -z "$(type csg_call)" ] && die "${0##*/}: csg_call not found"
tags=$(csg_call -l | awk -v name="$script" '($3==name){print $1,$2;exit}') || die "could not get tags"
helpmsg="$(csg_call $tags --help)" || die "${0##*/}: csg_call $tags --help failed"

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

echo
echo
echo
#do not parse xml options of functions files
[ -z "${script##function_*}" ] && exit 0

content="$(csg_call --cat $tags)" || die "${0##*/}: csg_call --cat $tags failed"
helpmsg="$(echo -e "$content" | sed -n -e '/^USES/d' -e '/csg_get_\(interaction_\)\?property/p')" || die "${0##*/}: sed failed"
echo "Used xml options:"
#set -x
#shell sciprts
if [ -n "$helpmsg" ] && [ -z "${script%%*.sh}" ]; then
  #1. trick to manage multiple csg_get_property per line by adding \n in the beginning
  #2. get value and their defaults and add link()() for txt2tags
  #3. rm quotes
  #4. adding itemize and link()() for txt2tags to create links later
  echo -e "$helpmsg" | \
  perl -n \
    -e 'if (/\(csg_get_(interaction_)?property\s+--allow-empty\s+(\S+?)\s*\)/) { print "$2 (default: empty)\n"; }
        elsif (/\(csg_get_(interaction_)?property\s+(\S+?)\s+(\S+?)\s*\)/) { print "$2 (default: $3)\n";}
        elsif (/\(csg_get_(interaction_)?property\s+(\S+?)\s*\)/) { print "$2\n";}
 	elsif (/csg_get_(interaction_)?property/) {die "Oho, I do NOT understand the line $_\n";}' | \
  sed -e 's/"//g' -e "s/'//g" | \
  perl -pe 's/^(\S+)/- link(PREFIX$1)($1)/;' -e 's/PREFIX([^c][^g])/cg.interaction.$1/;' -e 's/PREFIX//;'
fi
