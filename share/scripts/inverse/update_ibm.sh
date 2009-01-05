#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function update
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/} ATOM1 ATOM2
   echo Needs:  run_or_exit, \$source_wrapper, update_POT.pl
   exit 0
fi

if [ -z "$2" ]; then
   echo Missing argument >/dev/stderr
   exit 1
fi

if [ -z "$source_wrapper" ]; then
   echo source_wrapper not defined > /dev/stderr
   exit 1
fi
update_POT="$($source_wrapper --direct update_POT.pl)" || exit 1

if [ $(type -t run_or_exit) != "function" ]; then
   echo Could not find function run_or_exit > /dev/stderr
   exit 1
fi

local atom1=$1
local atom2=$2
run_or_exit --log log_update_POT_${atom1}_${atom2} $update_POT rdf_${atom1}_${atom2}_aim.xvg rdf_${atom1}_${atom2}.xvg delta_pot_${atom1}_${atom2}.d