#!/bin/bash 

die () {
  echo -e "$*" >&2
  exit 1
}

if [ "$1" = "--help" ]; then
  echo Usage: ${0##*/} xmlname
  echo Will convert xml to txt2tags file
  exit 0
fi

[ -z "$1" ] && die "${0##*/}: Missing argument"
[ -z "$(type -p csg_property)" ] && die "${0##*/}: csg_property not found"

xmlfile="$1"
[ -z "${CSGSHARE}" ] && die "${0##*/}: CSGSHARE not defined"
[ -f "${CSGSHARE}/xml/$xmlfile" ] || die "${0##*/}: Error, did not find ${CSGSHARE}/xml/$xmlfile"

#header lines
echo $xmlfile 
echo ${0##*/}
date
echo '%!includeconf: config.t2t'
echo

items="$(csg_property --file ${CSGSHARE}/xml/$xmlfile --path tags.item --print name --short)" || die "parsing xml failed"
for name in ${items}; do
  spaces="$(echo "${name//[^.]}" | sed -e 's/\./  /g')"
  echo "${spaces}- anchor(${name//\$})(**${name##*.}**)"
  desc="$(csg_property --file ${CSGSHARE}/xml/$xmlfile --path tags.item --filter "name=$name" --print desc --short)" || die "${0##*/}: Could not get desc for $name"
  echo ${desc}
done


