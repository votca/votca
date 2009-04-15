#! /bin/bash

if [ "$1" = "--help" ]; then
  echo This script implemtents the function initialize
  echo for the Inverse Boltzmann Method
  echo Usage: ${0##*/}
  echo Needs: SOURCE_WRAPPER, run_or_exit  
  exit 0
fi

RDF_to_POT="$($SOURCE_WRAPPER --direct RDF_to_POT.pl)" || exit 1
name=$($csg_get name)

if [ -f ../${name}.pot.in ]; then
  log Using given table for ${name}
  run_or_exit cp ../table_${name}.pot.in ${name}.pot.new 
else
  # RDF_to_POT.pl just does log g(r) + extrapolation
  log Using intial guess from RDF for ${name}
  run_or_exit $RDF_to_POT ${name}.dist.tgt ${name}.pot.new
fi

