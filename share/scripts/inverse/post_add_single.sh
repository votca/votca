#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script make all the post update with backup for single pairs 
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$SOURCE_WRAPPER
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

name=$($csg_get name)
tasklist=$($csg_get post_add) 
i=1
for task in $tasklist; do
  log "Doing $task for ${name}"
  mv ${name}.pot.new ${name}.pot.cur || die "${0##*/}: mv failed"
  cp ${name}.pot.cur ${name}.pot.${i} || die "${0##*/}: cp failed"
  script=$($SOURCE_WRAPPER postadd $task) || die "${0##*/}: $SOURCE_WRAPPER postadd $task failed"
  csg_get="$csg_get" bondtype="$bondbype" $script $1 || die "${0##*/}: $script $1 failed" 
  ((i++))
done
