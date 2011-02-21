#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the pressure update
   echo Usage: ${0##*/}
   exit 0
fi

for_all "non-bonded" \
  critical do_external table add '$(csg_get_interaction_property name).pot.cur $(csg_get_interaction_property name).dpot.new $(csg_get_interaction_property name).pot.new' 
