#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script make all the post update with backup for single pairs 
   echo Usage: ${0##*/} step_nr
   echo USES:  csg_get_interaction_property log mv die cp do_external run_or_exit
   echo NEEDS: name post_update
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

name=$(csg_get_interaction_property name)
tasklist=$(csg_get_interaction_property --allow-empty post_update) 
i=1
for task in $tasklist; do
  log "Doing $task for ${name}"
  run_or_exit mv ${name}.dpot.new ${name}.dpot.cur
  run_or_exit cp ${name}.dpot.cur ${name}.dpot.${i}
  run_or_exit do_external postupd "$task" "$1"
  ((i++))
done
