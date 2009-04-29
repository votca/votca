#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script make all the post update with backup 
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$source_wrapper, add_POT.pl
   exit 0
fi

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

update_single=$($SOURCE_WRAPPER --direct post_update_single.sh) || die "${0##*/}: $SOURCE_WRAPPER --direct post_update_single.sh"
for_all non-bonded ${update_single} $1

