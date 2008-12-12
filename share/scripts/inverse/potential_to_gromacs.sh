#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This is a wrapper to convert potential to gromacs
   echo Usage: ${0##*/} potential.d
   echo Needs: run_or_exit, \$scriptdir, table_to_xvg.pl
   exit 0
fi

if [ -z "$1" ]; then
   echo Give an argument > /dev/stderr
   exit 1
fi

if [ -z "$scriptdir" ]; then
   echo scriptdir not defined > /dev/stderr
   exit 1
fi

if [ ! -r "$scriptdir/RDF_to_POT.pl" ]; then
   echo Could not find RDF_to_POT.pl > /dev/stderr
   exit 1
fi

if [ $(type -t run_or_exit) != "function" ]; then
   echo Could not find function run_or_exit > /dev/stdout
   exit 1
fi

input="${1%.d}.d" 
output="${1%.d}.xvg" 
echo Convert $input to $output
run_or_exit --log log_table_to_xvg_${atom1}_${atom2} ${scriptdir}/table_to_xvg.pl $input $output || exit 1

