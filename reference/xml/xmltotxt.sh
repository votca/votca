#!/bin/bash 

unset -f die
die () {
  [[ -z "$CSGLOG" ]] || log --no-warn "$*"
  echo -e "$*" 1>&2
  log --no-warn "killing all processes...."
  #send kill signal to all process within the process groups
  kill 0
  exit 1
}

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

if [ ! -f "$1" ]; then
  echo "Inputfile '$1' not found" >&2
  exit 1
fi

#header lines
echo $1 > $2
date -r $1  +%F >> $2
echo >> $2
echo '%!includeconf: config.t2t' >> $2

for name in $(csg_property --file $1 --path tags.item --print name --short || die "parsing xml failed" ); do
  echo ": $name anchor(${name//\$})" >> $2
  echo "$(csg_property --file $1 --path tags.item --filter "name=$name" --print desc --short | sed -e 's/^[[:space:]]*//' -e '/^$/d' || die \"parsing xml failed\" )" >> $2
done


