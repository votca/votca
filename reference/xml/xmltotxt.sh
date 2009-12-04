#!/bin/bash 

die () {
  echo -e "$*" >&2
  exit 1
}

if [ "$1" = "--help" ]; then
  echo Usage: ${0##*/} infile 
  echo Will convert xml to txt2tags file
  exit 0
fi

[ -z "$1" ] && die "${0##*/}: Missing argument"
[ -z "$(type -p csg_property)" ] && die "${0##*/}: csg_property not found"
[ -f "$1" ] || die "Inputfile '$1' not found"

#header lines
echo $1 
echo ${0##*/}
date -r $1  +%F 
echo '%!includeconf: config.t2t'
echo

for name in $(csg_property --file $1 --path tags.item --print name --short || die "parsing xml failed" ); do
  echo ": $name anchor(${name//\$})"
  echo "$(csg_property --file $1 --path tags.item --filter "name=$name" --print desc --short | sed -e 's/^[[:space:]]*//' -e '/^$/d' || die \"parsing xml failed\" )"
done


