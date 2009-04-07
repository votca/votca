#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the pressure update
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$source_wrapper, add_POT.pl
   exit 0
fi

add_POT=$($SOURCE_WRAPPER --direct add_POT.pl) || die "${0##*/}: $SOURCE_WRAPPER --direct add_POT.pl failed" 

for_all non-bonded \
  ${add_POT} '$($csg_get type1)_$($csg_get type2).pot.cur $($csg_get type1)_$($csg_get type2).dpot.new $($csg_get type1)_$($csg_get type2).pot.new' >> $CSGLOG 2>&1
