#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function update
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/} ATOM1 ATOM2
   echo Needs:  run_or_exit, \$scriptdir, update_POT.pl
   exit 0
fi

if [ -z "$2" ]; then
   echo Missing argument >/dev/stderr
   exit 1
fi

if [ -z "$scriptdir" ]; then
   echo scriptdir not defined > /dev/stderr
   exit 1
fi

if [ ! -r "$scriptdir/update_POT.pl" ]; then
   echo Could not find RDF_to_POT.pl > /dev/stderr
   exit 1
fi

if [ $(type -t run_or_exit) != "function" ]; then
   echo Could not find function run_or_exit > /dev/stderr
   exit 1
fi

local atom1=$1
local atom2=$2
run_or_exit --log log_update_POT_${atom1}_${atom2} ../update_POT.pl rdf_${atom1}_${atom2}_aim.xvg rdf_${atom1}_${atom2}.xvg delta_pot_${atom1}_${atom2}.d