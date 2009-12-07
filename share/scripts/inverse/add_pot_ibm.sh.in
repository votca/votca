#! /bin/bash

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This script implemtents the pressure update

Usage: ${0##*/} step_nr

USES: do_external for_all run_or_exit csg_get_interaction_property

NEEDS: name
EOF
   exit 0
fi

check_deps "$0"

for_all "non-bonded" \
  run_or_exit do_external table add '$(csg_get_interaction_property name).pot.cur $(csg_get_interaction_property name).dpot.new $(csg_get_interaction_property name).pot.new' 
