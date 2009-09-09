#!/bin/bash 

if [ "$1" = "--help" ]; then
  echo Usage: ${0##*/} infile outfile
  echo Will convert xml to txt2tags file
  exit 0
fi

if [ -z "$2" ]; then
  echo "Missing argument" >&2
  exit 1
fi

if [ -z "$(type -p csg_property)" ]; then
  echo "csg_property not found" >&2
  exit 1
fi

#header lines
echo $1 > $2
date -r cgoptions.xml  +%F >> $2
echo >> $2

for name in $(csg_property --file $1 --path tags.item --print name --short ); do
  echo ": $name" >> $2
  echo "$(csg_property --file $1 --path tags.item --filter "name=$name" --print desc --short | sed -e 's/^[[:space:]]*//' -e '/^$/d' )" >> $2
done


