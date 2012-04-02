#!/bin/bash 

die () {
  echo -e "$*" >&2
  exit 1
}

add_heads() {
  local i head pattern again="no"
  for i in $items; do
    head=${i%.*}
    #for main node
    [ "$head" = "$i" ] && continue
    pattern=${head//$/\\$}
    #if head node is note there
    if [ -z "$(echo -e "$items" | grep -Ee "$pattern([[:space:]]+|\$)")" ]; then
      items="$(echo -e "$items\n$head")"
      again="yes"
      #echo "added $head from $i xx $pattern" >&2
      #break here to avoid double head nodes
      break
    fi
  done
  #we have to do that because items were changed
  [ "$again" = "yes" ] && add_heads
}

cut_heads() {
  local i new 
  [ -z "$1" ] && die "cut_heads: Missing argument"
  item="$1"
  hspace=0
  spaces=""
  for ((i=0;i<7;i++)); do
    new="${item#*.}"
    [ "$new" = "$item" ] && break
    item="$new"
    spaces+="  "
    let "hspace+=10"
  done
}

if [ "$1" = "--help" ]; then
  echo Usage: ${0##*/} xmlname
  echo Will convert xml to txt2tags file
  exit 0
fi

[ -z "$1" ] && die "${0##*/}: Missing argument"
[ -z "$(type -p csg_property)" ] && die "${0##*/}: csg_property not found"
[ -z "$(type -p csg_call)" ] && die "${0##*/}: csg_call not found"

xmlfile="$1"

# TODO: this is hacked
CSGSHARE="$(csg_call --show-share)/ctp"
[ -f "${CSGSHARE}/xml/$xmlfile" ] || die "${0##*/}: Error, did not find ${CSGSHARE}/xml/$xmlfile"

trunc=""
[ "$xmlfile" = "segments.xml" ] && trunc="segments."
[ "$xmlfile" = "mapping.xml" ] && trunc="mapping."
[ "$xmlfile" = "options.xml" ] && trunc="options."

#header lines
echo $xmlfile 
echo ${0##*/}
date
echo '%!includeconf: config.t2t'
echo

#get all items
items="$(csg_property --file ${CSGSHARE}/xml/$xmlfile --path tags.item --print name --short)" || die "parsing xml failed"
#check if a head node is missing
#add_heads
#sort them
#items="$(echo -e "$items" | sort -u)"
#echo "$items"

for name in ${items}; do
  #cut the first 3 heads of the item
  cut_heads "$name"
  head="$(echo $name | sed -e 's/'${item}'//')"
  #echo $head $name $item
  desc="$(csg_property --file ${CSGSHARE}/xml/$xmlfile --path tags.item --filter "name=$name" --print desc --short)" || die "${0##*/}: Could not get desc for $name"

  if [ "$name" = "sectionlink" ]; then
     link="$name($desc)"
  else
     echo " | hspace($hspace) target(${trunc}${name})(${item})  | ${desc}"
  fi
done

if  [ ! -z "$link" ]; then 
  echo  "Return to the descritpion of $link." 
fi

