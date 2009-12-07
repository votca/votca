#! /bin/bash

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This script make all the post update with backup

Usage: ${0##*/} step_nr

USES:  do_external die for_all

NEEDS:
EOF
   exit 0
fi

check_deps "$0"

[[ -n "$1" ]] || die "${0##*/}: Missing argument"

for_all "non-bonded" do_external post update_single ${1}

