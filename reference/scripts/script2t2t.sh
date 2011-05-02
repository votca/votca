#! /bin/bash

die() {
  echo -e "$*" >&2
  exit 1
}

assert() {
  local x pipestatus=${PIPESTATUS[*]}
  for x in $pipestatus ; do
      [[ $x -eq 0 ]] || die "$@"
  done
}

[ -z "$1" ] && die "${0##*/}: Missing argument"
script="$1"
shift

[ -z "$(type csg_call)" ] && die "${0##*/}: csg_call not found"
tags=$(csg_call -l | awk -v name="$script" '($3==name){print $1,$2;exit}') 
assert "could not get tags"
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
#-remove script version line
#-space at begining and end are removed
#-optionlines (-[^ ]) -> - ``--option`` text
#-usageline -> usage: ``code``
#-examplelines (^*) -> - ``line``
#-extra empty line before new section to close txt2tags itemize (2 empty lines)
#-link external packages
echo -e "$helpmsg" | sed \
   -e '1,/^please submit/d' \
   -e "/^$script/d" \
   -e 's/^[[:space:]]*//' \
   -e 's/[[:space:]]*$//' \
   -e '/^-[^ ]/s/ \{2\}/`` /' \
   -e '/^-[^ ].*``/s/^/- ``/' \
   -e 's/^\(Usage:[[:space:]]*\)\(.*\)$/\1``\2``/' \
   -e '/^\* /s/\( \{2\}\|$\)/`` /' \
   -e '/^\*.*``/s/^\*[[:space:]]*/- ``/' \
   -e 's/^\(Examples\|Usage\):/\n&/' \
   -e 's/^\(Allowed\|Trajectory\|Specific\) options:/\n&/' \
   -e 's/^\(Used external packages:\)[[:space:]]*\(.*\)$/\1 ref(progpack.\2)(\2)/'

assert "sed 1 failed"
#do not parse xml options of functions files
[ -z "${script##functions_*}" ] && exit 0

content="$(csg_call --cat $tags)" || die "${0##*/}: csg_call --cat $tags failed"
helpmsg="$(echo -e "$content" | sed -n -e '/csg_get_\(interaction_\)\?property/p')" 
assert "${0##*/}: sed 2 failed"

#no properties found
[ -n "$helpmsg" ] || exit 0

echo
echo "Used xml options:"
#set -x

#shell sciprts
#1. trick to manage multiple csg_get_property per line by adding \n in the beginning
#2. get value and possibly their defaults
#3. rm quotes
#4. adding itemize and link()() (see txt2tags config.t2t) and add cg.interaction to the link target for interaction properties
#5. remove duplicated
echo -e "$helpmsg" | \
sed 's/csg_get_/\n&/g' | \
if [ -z "${script%%*.sh}" ]; then
  perl -n \
    -e 'if (/csg_get_(interaction_)?property\s+--allow-empty\s+(\S+?)\s*\)/) { print "$2 (default: empty)\n"; }
        elsif (/csg_get_(interaction_)?property\s+(\S+?)\s+(\S+?)\s*\)/) { print "$2 (default: $3)\n";}
        elsif (/csg_get_(interaction_)?property\s+(\S+?)\s*\)/) { print "$2\n";}
 	elsif (/csg_get_(interaction_)?property/) {die "Oho, I do NOT understand the line $_\n";}'
elif [ -z "${script%%*.pl}" ]; then
  perl -n \
    -e 'if (/csg_get_(interaction_)?property\s*\(\s*(\S+?)\s*,\s*(\S+?)\s*\)/) { print "$2 (default: $3)\n";}
        elsif (/csg_get_(interaction_)?property\s*\(\s*(\S+?)\s*\)/) { print "$2\n";}
 	elsif (/csg_get_(interaction_)?property/) {die "Oho, I do NOT understand the line $_\n";}'
else
  die "Don't know how to handle script ${script}"
fi | \
sed -e 's/"//g' -e "s/'//g" | \
perl -pe 's/^(\S+)/- link(PREFIX$1)($1)/;' -e 's/PREFIX([^c][^g])/cg.non-bonded.$1/;' -e 's/PREFIX//;' | \
sort -u

assert "${0##*/}: sed 3 failed"
