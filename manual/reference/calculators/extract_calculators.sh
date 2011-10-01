#!/bin/bash -e


VOTCASHARE="$(csg_call --show-share)"
calculators="$(ctp_run --list | sed -ne 's/^\s\+\([a-z]*\)\s*\(.*\)/\1/p')"
texfile="$PWD/calculators.tex"

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
  for ((i=0;i<5;i++)); do
    new="${item#*.}"
    [ "$new" = "$item" ] && break
    item="$new"
    spaces+="  "
    let "hspace+=10"
  done
}

# get the loop for all calculators
echo calculators: $calculators
rm -f $texfile

for calculator in ${calculators}; do
  xmlfile=${VOTCASHARE}/ctp/xml/$calculator.xml
  calculator_description="$(csg_property --file $xmlfile --path $calculator --print description --short)"
  calculator_sectionlabel="$(csg_property --file $xmlfile --path $calculator --print sectionlabel --short)"

  echo "\subsection{$calculator}" >> $texfile
  echo "\label{calc:$calculator}" >> $texfile
  echo "%" >> $texfile
  echo $calculator_description >> $texfile
  echo "%" >> $texfile
  echo "\rowcolors{1}{invisiblegray}{white}" >> $texfile
  echo "\begin{longtable}{m{3cm}|m{11cm}}" >> $texfile

  # including XML files from the reference section.
  #echo "\input{$calculator.xml}" >> $texfile

  echo; echo " "$calculator [sec:$calculator_sectionlabel]: $calculator_description

  items="$(csg_property --file $xmlfile --path $calculator.item --print name --short)" || die "parsing xml failed"
  echo "   "properties: $items

  for name in ${items}; do  
    cut_heads "$name"
    head="$(echo $name | sed -e 's/'${item}'//')"
    description="$(csg_property --file $xmlfile --path $calculator.item --filter "name=$name" --print description --short)" || die "${0##*/}: Could not get desc for $name"
    echo "      "$head $name $item $description
    echo " \hspace{${hspace}pt} \hypertarget{$calculator.${trunc}${name}}{${item}}  & ${description} \\\\" >> $texfile
  done

  echo "\end{longtable}" >> $texfile
  echo "%" >> $texfile 
  if [ -n "$calculator_sectionlabel" ]; then
    echo "Return to the descritpion of \slink{$calculator_sectionlabel}{\texttt{$calculator}}" >> $texfile
  fi
  echo >> $texfile

done




