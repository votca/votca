#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script make all the post update with backup for single pairs 
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$source_wrapper, add_POT.pl
   exit 0
fi

add_POT=$($SOURCE_WRAPPER --direct add_POT.pl) || die "${0##*/}: $SOURCE_WRAPPER --direct add_POT.pl"

type1=$($csg_get type1)
type2=$($csg_get type2)
tasklist=$($csg_get post_add) 
#we type1_type2.dpot.new
i=1
for task in $tasklist; do
  log "Doing $task for $type1 $type2"
  mv ${type1}_${type2}.dpot.new ${type1}_${type2}.dpot.cur || die "${0##*/}: mv failed"
  cp ${type1}_${type2}.dpot.cur ${type1}_${type2}.dpot.${i} || die "${0##*/}: cp failed"
  script=$($SOURCE_WRAPPER postadd $task) || die "${0##*/}: $SOURCE_WRAPPER postadd $task failed"
  csg_get="$csg_get" bondtype="$bondbype" "$script" >> $CSLOG 2>&1 || die "${0##*/}: $script failed" 
  ((i++))
done
