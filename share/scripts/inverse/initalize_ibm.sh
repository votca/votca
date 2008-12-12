#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function initialize
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo Needs: \$atoms, run_or_exit, \$scriptdir, RDF_to_POT.pl
   exit 0
fi

if [ -z "$scriptdir" ]; then
   echo scriptdir not defined > /dev/stderr
   exit 1
fi

if [ ! -r "$scriptdir/RDF_to_POT.pl" ]; then
   echo Could not find RDF_to_POT.pl > /dev/stderr
   exit 1
fi

if [ ${#atoms[@]} -eq 0 ]; then
   echo Could not get atomname for \$atoms > /dev/stderr
   exit 1
fi

if [ $(type -t run_or_exit) != "function" ]; then
   echo Could not find function run_or_exit > /dev/stderr
   exit 1
fi

for ((i=0;i<${#atoms[*]};i++)); do
   atom1=${atoms[$i]}
   for ((j=$i;j<${#atoms[*]};j++)); do
      atom2=${atoms[$j]}
      if [ -f ../table_${atom1}_${atom2}_guess.d ]; then
         echo Using given table for $atom1-$atom2
         cp ../table_${atom1}_${atom2}_guess.d table_${atom1}_${atom2}_new.d || exit 1
      else
         # RDF_to_POT.pl just does log g(r) + extrapolation
         echo Using intial guess from RDF for ${atom1}-${atom2}
         run_or_exit  --log log_RDF_to_POT_${atom1}_${atom2} $scriptdir/RDF_to_POT.pl rdf_${atom1}_${atom2}_aim.xvg table_${atom1}_${atom2}_new.d || exit 1
      fi
   done
done
