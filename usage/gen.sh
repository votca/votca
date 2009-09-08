#!/bin/bash
xml=$1

echo "\\begin{itemize}" 

for name in $(csg_property --file $xml --path tags.item --print name --short | sort); do
  echo "\\item \\textbf{$name} \\\\" 
  echo "$(csg_property --file $xml --path tags.item --filter \"name=$name\" --print desc --short | awk '{if(NF>0 && !first) print $0; else first=true;}' )"
done

echo "\\end{itemize}"

