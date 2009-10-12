#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script make all the post update with backup 
   echo Usage: ${0##*/} step_nr
   echo USES:  do_external die for_all
   echo NEEDS:
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

for_all "non-bonded" do_external post update_single ${1}

