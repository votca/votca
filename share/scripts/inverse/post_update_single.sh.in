#! /bin/bash

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This script make all the post update with backup for single pairs

Usage: ${0##*/} step_nr

USES:  csg_get_interaction_property log mv die cp do_external run_or_exit

NEEDS: name inverse.post_update
EOF
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

name=$(csg_get_interaction_property name)
tasklist=$(csg_get_interaction_property --allow-empty inverse.post_update) 
i=1
for task in $tasklist; do
  log "Doing $task for ${name}"
  run_or_exit mv ${name}.dpot.new ${name}.dpot.cur
  run_or_exit cp ${name}.dpot.cur ${name}.dpot.${i}
  do_external postupd "$task" "$1"
  ((i++))
done
