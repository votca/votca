#!/bin/bash
# Script to generate .fig dependency file
# 

usage() {
echo $@ 1>&2
exit 1
}

# check if there is only one input parameter
[ $# -ne 1 ] && usage "Error: one parameter is expected - name of the .fig file"

FILENAME="$1"
BASENAME=${FILENAME%.fig}

# check if argument is actually a .fig file
if [ "X${FILENAME}" = "X${BASENAME}" ]
then
   usage "Error: .fig file must be supplied"
fi

DEPS=$(sed -n '/.*eps\|png/ s/[[:blank:]]*[0-9][[:blank:]]*//p' $FILENAME | \
   awk -v ORS=' ' '{print $0}')

echo ${BASENAME}.eps : $DEPS 
