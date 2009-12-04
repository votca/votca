#! /bin/bash

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version @version@
This is a intialize stuff for gromacs

Usage: ${0##*/}

USES: cp die

NEEDS:
EOF
   exit 0
fi

check_deps "$0"

cp ../conf.gro confout.gro || die "cp ../conf.gro confout.gro"
