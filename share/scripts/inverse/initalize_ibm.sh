#! /bin/bash

if [ "$1" = "--help" ]; then
  echo This script implemtents the function initialize
  echo for the Inverse Boltzmann Method
  echo Usage: ${0##*/}
  echo Needs: SOURCE_WRAPPER, run_or_exit  
  exit 0
fi

if [ -z "${SOURCE_WRAPPER}" ]; then
  echo SOURCE_WRAPPER not defined > /dev/stderr
  exit 1
fi

if [ "$(type -t run_or_exit)" != "function" ]; then
  echo Could not find function run_or_exit > /dev/stderr
  exit 1
fi

RDF_to_POT="$($SOURCE_WRAPPER --direct RDF_to_POT.pl)" || exit 1

if [ -f ../table_${type1}_${type2}_guess.d ]; then
   echo Using given table for $type1-$type2
   cp ../table_${type1}_${type2}_guess.d table_${type1}_${type2}_new.d || exit 1
else
   # RDF_to_POT.pl just does log g(r) + extrapolation
   echo Using intial guess from RDF for ${type1}-${type2}
   run_or_exit  --log log_RDF_to_POT_${type1}_${type2} $RDF_to_POT rdf_${type1}_${type2}_aim.xvg table_${type1}_${type2}_new.d || exit 1
fi

