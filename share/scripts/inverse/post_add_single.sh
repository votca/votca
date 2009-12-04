#! /bin/bash

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This script make all the post update with backup for single pairs 

Usage: ${0##*/} step_nr

USES:  csg_get_interaction_property log mv cp do_external run_or_exit die

NEEDS: name inverse.post_add
EOF
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

name=$(csg_get_interaction_property name)
tasklist=$(csg_get_interaction_property --allow-empty inverse.post_add) 
i=1
for task in $tasklist; do
  log "Doing $task for ${name}"
  run_or_exit mv ${name}.pot.new ${name}.pot.cur
  run_or_exit cp ${name}.pot.cur ${name}.pot.${i}
  do_external postadd "$task" "$1"
  ((i++))
done
