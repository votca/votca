#! /bin/bash

if [ "$1" = "--help" ]; then
  cat <<EOF
Usage: ${0##*/} settings.xml

This script needs \$CSGSHARE to find inverse.sh,
so export it before or give a file where it is defined

Examples:
         CSGSHARE=\$HOME/csg/csg/share ${0} xml file
         ${0##*/} xml file

EOF
  exit 0
fi

if [ -z "$1" ]; then
  echo Error: no xml file given 1>&2
  echo Try ${0##*/} --help for help 1>&2
  exit 1
fi

#if not defined try to parse agrument
if [ -z "${CSGSHARE}" ]; then
   echo Error: Environment values CSGSHARE not defined 1>&2
   echo Try ${0##*/} --help for help 1>&2
   exit 1
fi

echo \$CSGSHARE is ${CSGSHARE}

if [ -f ${CSGSHARE}/scripts/inverse/inverse.sh ]; then
   exec ${CSGSHARE}/scripts/inverse/inverse.sh $1
   exit 0
else
   echo ${0##*/}: Could not run ${CSGSHARE}/scripts/inverse/inverse.sh $1
   exit 1
fi
