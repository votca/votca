#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function initialize
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo Needs: \$atoms, run_or_exit, $\source_wrapper.sh
   exit 0
fi

if [ -z "$source_wrapper" ]; then
   echo source_wrapper not defined > /dev/stderr
   exit 1
fi
RDF_to_POT="$($source_wrapper --direct RDF_to_POT.pl)" || exit 1

if [ ${#atoms[@]} -eq 0 ]; then
   echo Could not get atomname for \$atoms > /dev/stderr
   exit 1
fi

if [ $(type -t run_or_exit) != "function" ]; then
   echo Could not find function run_or_exit > /dev/stderr
   exit 1
fi

for ((a1=0;a1<${#atoms[*]};a1++)); do
   atom1=${atoms[$a1]}
   for ((a2=$a1;a2<${#atoms[*]};a2++)); do
      atom2=${atoms[$a2]}
      if [ -f ../table_${atom1}_${atom2}_guess.d ]; then
         echo Using given table for $atom1-$atom2
         cp ../table_${atom1}_${atom2}_guess.d table_${atom1}_${atom2}_new.d || exit 1
      else
         # RDF_to_POT.pl just does log g(r) + extrapolation
         echo Using intial guess from RDF for ${atom1}-${atom2}
         run_or_exit  --log log_RDF_to_POT_${atom1}_${atom2} $RDF_to_POT rdf_${atom1}_${atom2}_aim.xvg table_${atom1}_${atom2}_new.d || exit 1
      fi
   done
done
