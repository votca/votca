#! /bin/bash

if [ "$1" = "--help" ]; then
   echo "This script implemtents smoothing of the potential update (.dpot)"
   echo "Usage: ${0##*/} step_nr"
   echo "Needs:  run_or_exit, \$SOURCE_WRAPPER"
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

name=$($csg_get name)
tmpfile=$(mktemp ${name}.XXX) || die "mktemp failed"
iterations=$($csg_get post_update_options.smooth.iterations)  

smooth=$($SOURCE_WRAPPER table smooth) || die "${0##*/}: $SOURCE_WRAPPER table smooth failed"

cp ${name}.dpot.cur $tmpfile
log "doing $iterations smoothing iterations"

for((i=0;i<$iterations;i++)); do
  run_or_exit ${smooth} $tmpfile ${name}.dpot.new
  cp ${name}.dpot.new $tmpfile
done

