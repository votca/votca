#! /bin/bash

if [ "$1" = "--help" ]; then
  echo This script implemtents the function initialize
  echo for the Inverse Boltzmann Method
  echo Usage: ${0##*/}
  echo Needs: \$SOURCE_WRAPPER, run_or_exit  
  exit 0
fi

RDF_to_POT="$($SOURCE_WRAPPER --direct RDF_to_POT.pl)" || exit 1
name=$($csg_get name)

if [ -f ../${name}.pot.in ]; then
  log "Using given table ${name}.pot.in for ${name}"
  min=$($csg_get min )
  max=$($csg_get max )
  step=$($csg_get step )
  run_or_exit csg_resample --in ../${name}.pot.in --out ${name}.pot.new --grid ${min}:${step}:${max}
else
  # RDF_to_POT.pl just does log g(r) + extrapolation
  log "Using intial guess from RDF for ${name}"
  run_or_exit ${RDF_to_POT} ${name}.dist.tgt ${name}.pot.new
fi

