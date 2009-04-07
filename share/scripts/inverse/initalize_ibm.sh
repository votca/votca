#! /bin/bash

if [ "$1" = "--help" ]; then
  echo This script implemtents the function initialize
  echo for the Inverse Boltzmann Method
  echo Usage: ${0##*/}
  echo Needs: SOURCE_WRAPPER, run_or_exit  
  exit 0
fi

RDF_to_POT="$($SOURCE_WRAPPER --direct RDF_to_POT.pl)" || exit 1
type1=$($csg_get type1)
type2=$($csg_get type2)

if [ -f ../${type1}_${type2}.pot.in ]; then
  log Using given table for $type1-$type2
  run_or_exit cp ../table_${type1}_${type2}pot.in ${type1}_${type2}.pot.new 
else
  # RDF_to_POT.pl just does log g(r) + extrapolation
  log Using intial guess from RDF for ${type1}-${type2}
  run_or_exit $RDF_to_POT ${type1}_${type2}.dist.tgt ${type1}_${type2}.pot.new
fi

