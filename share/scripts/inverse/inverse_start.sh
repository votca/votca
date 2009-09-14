#! /bin/bash

#NEEDS: cg.inverse.scriptdir cg.inverse.log_file cg.inverse.restart_file
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

#find source_wrapper.pl
SOURCE_WRAPPER="${CSGINVERSE}/source_wrapper.pl"
[[ -x "${SOURCE_WRAPPER}" ]] || die "Could not find source_wrapper.pl"
export SOURCE_WRAPPER

function_file=$($SOURCE_WRAPPER functions common) || die "$SOURCE_WRAPPER functions common failed"
#die() is overwritten here
source ${function_file} || exit 1
unset function_file

CSGSCRIPTDIR="$(csg_get_property cg.inverse.scriptdir)" 
#scriptdir maybe contains $PWD or something
eval CSGSCRIPTDIR=$CSGSCRIPTDIR
[[ -d "$CSGSCRIPTDIR" ]] || die "CSGSCRIPTDIR '$CSGSCRIPTDIR' is not a dir"
export CSGSCRIPTDIR

CSGLOG="$(csg_get_property cg.inverse.log_file)"
CSGLOG="$PWD/$CSGLOG"
export CSGLOG

CSGRESTART="$(csg_get_property cg.inverse.restart_file)"
export CSGRESTART

