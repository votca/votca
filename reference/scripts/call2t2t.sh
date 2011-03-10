#! /bin/bash

if ! which csg_call > /dev/null 2>&1; then
  echo csg_call not found
  exit 1
fi

echo "csg_table"
echo "make"
date
echo
echo '%!includeconf: config.t2t'
echo
echo "|| Key1 | Key2 | Scriptname"
csg_call -l | \
sed -e '1d' -e '/\.\(sh\|pl\)/!d' -e '/^functions/d' | \
awk '{printf "| %s | %s | ref(%s)(%s) |\n",$1,$2,$3,$3}'
