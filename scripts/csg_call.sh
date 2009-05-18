#!/bin/bash

#echo "$CSGSHARE/scripts/inverse/inverse_start.sh"
#source "$CSGSHARE/scripts/inverse/inverse_start.sh"  "$@" || exit 1

export CSGINVERSE=$CSGSHARE/scripts/inverse
export SOURCE_WRAPPER=$CSGINVERSE/source_wrapper.pl
export CSGLOG=log.tmp
source $($SOURCE_WRAPPER functions common) || die "$SOURCE_WRAPPER functions common failed" 

script=$($SOURCE_WRAPPER $1 $2)

echo running $script

do_external $*