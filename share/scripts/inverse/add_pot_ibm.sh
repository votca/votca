#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the pressure update
   echo Usage: ${0##*/} step_nr
   echo USES: \$SOURCE_WRAPPER die for_all run_or_exit csg_get_interaction_property
   echo NEEDS: name
   exit 0
fi

check_deps "$0"

add_POT=$($SOURCE_WRAPPER add pot) || die "${0##*/}: $SOURCE_WRAPPER add pot failed" 

for_all non-bonded \
  run_or_exit ${add_POT} '$(csg_get_interaction_property name).pot.cur $(csg_get_interaction_property name).dpot.new $(csg_get_interaction_property name).pot.new' 
