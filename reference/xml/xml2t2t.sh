#!/bin/bash 

die () {
  echo -e "$*" >&2
  exit 1
}

cut_heads() {
  local i new 
  [ -z "$1" ] && die "cut_heads: Missing argument"
  item="$1"
  spaces=""
  for ((i=0;i<3;i++)); do
    new="${item#*.}"
    [ "$new" = "$item" ] && break
    item="$new"
    spaces+="  "
  done
}

if [ "$1" = "--help" ]; then
  echo Usage: ${0##*/} xmlname
  echo Will convert xml to txt2tags file
  exit 0
fi

[ -z "$1" ] && die "${0##*/}: Missing argument"
[ -x "${CSG_PROPERTY}" ] || die "${0##*/}: variable CSG_PROPERTY (value ${CSG_PROPERTY}) is wrong"

xmlfile="$1"
[ -f "$xmlfile" ] || die "${0##*/}: Error, did not find $xmlfile"

#header lines
echo ${xmlfile##*/}
echo ${0##*/}
date
echo '%!includeconf: config.t2t'
echo
echo **Please mind that dots in xml tags have to replaced by subtags, e.g. x.y has to be converted to x with subtag y.**
echo

#get all items
#perl does this:
#-keep only xml tags
#-print the xml tags tree
#-remove leading dot
#-delete desc tags themself
items="$(perl -ne 'while (/<([^!].*?)>/g) {print "$1\n";}' ${xmlfile} | \
  perl -ne 'chomp;if (/^\/(.*)/) {print "$line\n";$line =~ s/\.$1$//;}else{$line.=".$_";};' | \
  sed -e 's/^\.//' -e '/\.DESC$/d')"

#sort them
items="$(echo -e "$items" | sort -u)"
#echo "$items"
for name in ${items}; do
  #cut the first 3 heads of the item, as maximum deepth of latex itemize is 3
  cut_heads "$name" #returns $item and $spaces
  desc="$(${CSG_PROPERTY} --file $xmlfile --path $name.DESC --print . --short)"
  echo -n "${spaces}- target(${name})(**${item}**) "
  default="$(${CSG_PROPERTY} --file $xmlfile --path $name --print . --short | \
    sed -e '/^[[:space:]]*$/d' -e 's/^[[:space:]]*//' -e 's/[[:space:]]$//')"
  [[ -n $default ]] && desc="$desc (default $default)"
  echo ${desc}
done


