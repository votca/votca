#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script make all the post update with backup 
   echo Usage: ${0##*/} step_nr
   echo USES:  \$SOURCE_WRAPPER die for_all
   echo NEEDS:
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

update_single=$($SOURCE_WRAPPER postadd single) || die "${0##*/}: $SOURCE_WRAPPER postadd single failed"
for_all non-bonded ${update_single} ${1}

