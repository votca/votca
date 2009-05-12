#! /bin/bash

#this script just contain all startup check for inverse.sh to 
#make inverse.sh more readable

#for now, we will replace this function later
die(){ echo "$*" >&2; exit 1; }

[[ -n "$1" ]] || die "Error: Missing xml file"

if [ -f "./$1" ]; then
  export CSGXMLFILE="${PWD}/${1}"
else
  die "Error: file could not read"
fi

#check for CSGSHARE 
[[ -n "${CSGSHARE}" ]] || die "Error: CSGSHARE not definded"
[[ -d "${CSGSHARE}" ]] || die "CSGSHARE '$CSGSHARE' is not a dir"

CSGINVERSE="${CSGSHARE}/scripts/inverse"
[[ -d "${CSGINVERSE}" ]] || die "CSGSHARE/scripts/inverse is not a dir"
export CSGINVERSE

#we need csg_property
[[ -n "$(type -p csg_property)" ]] || die "Error: csg_property not found, check your PATH"

CSGSCRIPTDIR="$(csg_property --file $CSGXMLFILE --path cg.inverse.scriptdir --short --print .)" || \
  die "csg_property --file $CSGXMLFILE --path cg.inverse.scriptdir --short --print . failed" 
#scriptdir maybe contains $PWD or something

eval CSGSCRIPTDIR=$CSGSCRIPTDIR
[[ -d "$CSGSCRIPTDIR" ]] || die "CSGSCRIPTDIR '$CSGSCRIPTDIR' is not a dir"
export CSGSCRIPTDIR

#find source_wrapper.pl
#first local scriptdir. then CSGINVERSE
if [ -f "${CSGSCRIPTDIR}/source_wrapper.pl" ]; then
   SOURCE_WRAPPER="${CSGSCRIPTDIR}/source_wrapper.pl"
elif [ -f "${CSGINVERSE}/source_wrapper.pl" ]; then
   SOURCE_WRAPPER="${CSGINVERSE}/source_wrapper.pl"
else
   die "Could not find source_wrapper.pl"
fi
export SOURCE_WRAPPER

CSGLOG="$(csg_property --file $CSGXMLFILE --path cg.inverse.log_file --short --print .)" || die "Could not get logfile"
[[ -n "${CSGLOG}" ]] || die "logfile is empty"
CSGLOG="$PWD/$CSGLOG"
export CSGLOG

function_file=$($SOURCE_WRAPPER functions common) || die "$SOURCE_WRAPPER functions common failed"
#die() is overwritten here
source ${function_file} || exit 1
unset function_file

