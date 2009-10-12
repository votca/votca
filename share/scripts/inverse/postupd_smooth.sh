#! /bin/bash

if [ "$1" = "--help" ]; then
   echo "This script implemtents smoothing of the potential update (.dpot)"
   echo "Usage: ${0##*/} step_nr"
   echo USES:  die csg_get_interaction_property mktemp do_external cp log run_or_exit
   echo NEEDS: name post_update_options.smooth.iterations
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

name=$(csg_get_interaction_property name)
tmpfile=$(mktemp ${name}.XXX) || die "mktemp failed"
iterations=$(csg_get_interaction_property post_update_options.smooth.iterations 1)  

run_or_exit cp ${name}.dpot.cur $tmpfile
log "doing $iterations smoothing iterations"

for((i=0;i<$iterations;i++)); do
  run_or_exit do_external table smooth $tmpfile ${name}.dpot.new
  run_or_exit cp ${name}.dpot.new $tmpfile
done

