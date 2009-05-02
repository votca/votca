#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script make all the post update with backup for single pairs 
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$source_wrapper, add_POT.pl
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

name=$($csg_get name)
tasklist=$($csg_get post_update) 
i=1
for task in $tasklist; do
  log "Doing $task for ${name}"
  mv ${name}.dpot.new ${name}.dpot.cur || die "${0##*/}: mv failed"
  cp ${name}.dpot.cur ${name}.dpot.${i} || die "${0##*/}: cp failed"
  script=$($SOURCE_WRAPPER postupd $task) || die "${0##*/}: $SOURCE_WRAPPER postupd $task failed"
  csg_get="$csg_get" bondtype="$bondbype" $script $1 || die "${0##*/}: $script $1 failed"
  ((i++))
done
