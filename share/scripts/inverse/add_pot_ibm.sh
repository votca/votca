#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the pressure update
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$source_wrapper, add_POT.pl
   exit 0
fi

add_POT=$($SOURCE_WRAPPER --direct add_POT.pl) || die "${0##*/}: $SOURCE_WRAPPER --direct add_POT.pl failed" 

for_all non-bonded \
  logrun ${add_POT} '$($csg_get name).pot.cur $($csg_get name).dpot.new $($csg_get name).pot.new' 
