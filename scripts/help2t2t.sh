#! /bin/bash -e

die() {
  echo -e "$*" >&2
  exit 1
}

[ -z "$1" ] && die "Usage: ${0##*/} prog"
prog="$1"

[ -x "./$prog" ] || die "${0##*/}: ./$prog not found"

#trash the header
helpmsg="$(./$prog --help | sed -e '1,10d')" || die "$prog --help failed"

#first block is descprition
desc="$(echo "$helpmsg" | sed -n '1,/^[[:space:]]*$/p')" || die "parse of desc failed"
[ -z "$desc" ] && die "Failed to fetch desc"
helpmsg="$(echo "$helpmsg" | sed '1,/^[[:space:]]*$/d')" || die "cut of help msg failed"

#second block can be usage line
usage="$(echo "$helpmsg" | sed -n '1s/Usage:[[:space:]]*\(.*\)$/\1/p')" || die "parse of usage failed"
if [ -z "$usage" ]; then
  usage="**$prog** \"\"[\"\"//OPTION//\"\"]\"\" \"\"[\"\"//OPTIONPARAMETERS//\"\"]\"\""
else
  usage="$(echo "$usage" | sed -e 's/^/**/' -e 's/ /** /' -e 's/\([][]\)/""\1""/g')" || \
    die "parse part 2 of usage failed"
  helpmsg="$(echo "$helpmsg" | sed '1,/^[[:space:]]*$/d')" || die "cut of help msg failed"
fi

#try to find examples block
exam="$(echo "$helpmsg" | sed -n '/^Examples:/,/^[[:space:]]*$/p')" || die "parse of exam failed"
if [ -n "$exam" ]; then
  exam="$(echo "$exam" | \
    sed -e '1d' \
      -e '/^\* /s/\( \{2\}\|$\)/``/' \
      -e '/^\*.*``/s/^\*[[:space:]]*/- ``/')" ||
    die "parse part 2 of exam failed"
  helpmsg="$(echo "$helpmsg" | sed '/^Examples:/,/^[[:space:]]*$/d')" || die "cut of help msg failed"
fi

#write t2t file
cat <<EOF
$prog
${0##*/}
$(date)

= NAME =

$prog - Part of the VOTCA package

= SYNOPSIS =

$usage

= DESCRIPTION =

$desc

Please visit the program site at __http://www.votca.org__.

= OPTIONS =
EOF
#Here comes the magic:
#-space at begining and end are removed
#-optionlines (-[^ ]) -> - ``--option`` text
#-usageline -> usage: ``code``
#-extra empty line before new section to close itemize
#-examplelines (^*) -> - ``line``
#-remove NEEDS and OPTINAL, ... line
#-added newline before new option block
echo -e "$helpmsg" | sed \
   -e 's/^[[:space:]]*//' \
   -e 's/[[:space:]]*$//' \
   -e '/^-[^ ]/s/ \{2,\}/**\n/' \
   -e '/^-[^ ].*\*\*/s/^/: **/' \
   -e '/^\* /s/\( \{2\}\|$\)/``/' \
   -e '/^\*.*``/s/^\*[[:space:]]*/- ``/' \
   -e '/^\(NEEDS\|OPTIONAL\|USES\|PROVIDES\)/d' \
   -e 's/^\(Allowed\|Trajectory\|Specific\|Mapping\) options:/\n&/'

if [ -n "$exam" ]; then
  cat <<EOF
= EXAMPLES =
$exam
EOF
fi
cat <<EOF


= AUTHORS =

Written and maintained by the VOTCA Development Team <devs@votca.org>

This Manual Page was extracted from '$prog --help', 
then converted to this format by [txt2tags http://txt2tags.org] !


= COPYRIGHT =

Copyright 2009-2010 The VOTCA Development Team

This is free software; see the source for copying conditions. There is
NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.
EOF
