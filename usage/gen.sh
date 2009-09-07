#!/bin/bash
file="mapping_tags.tex"
xml="mapping.xml"
echo "" > $file

echo "\\begin{itemize}" >> $file

for name in $(csg_property --file $xml --path tags.item --print name --short | sort); do
  echo "\\item \\textbf{$name} \\\\" >> $file
  echo "$(csg_property --file $xml --path tags.item --filter \"name=$name\" --print desc --short | awk '{if(NF>0 && !first) print $0; else first=true;}' )" >> $file
done

echo "\\end{itemize}" >> $file

