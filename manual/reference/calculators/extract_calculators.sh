#!/bin/bash -e

VOTCASHARE="$(csg_call --show-share)"
calculators="$(ctp_run --list | sed -ne 's/^\s\+\([a-z]*\)\s*\(.*\)/\1/p')"
texfile="$PWD/calculators.tex"

rm -f $texfile; touch $texfile

# loop over all calculators
for calculator in ${calculators}; do

  xmlfile=${VOTCASHARE}/ctp/xml/$calculator.xml

  if [ ! -f "$xmlfile" ]; then 
    continue
  fi

  votca_property --file $xmlfile --format TEX >> $texfile

done




