#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This is a wrapper to convert potential to gromacs
   echo Usage: ${0##*/} potential.d
   echo Needs: run_or_exit, \$source_wrapper, spline
   exit 0
fi

if [ -z "$1" ]; then
   echo Give an argument > /dev/stderr
   exit 1
fi

if [ -z "$source_wrapper" ]; then
   echo source_wrapper not defined > /dev/stderr
   exit 1
fi
table_to_xvg="$($source_wrapper --direct table_to_xvg.pl)" || exit 1

for exe in spline; do
   if [ -z $(type -p $exe) ]; then
      echo Could not find $exe > /dev/stderr
      exit 1
   fi
done

if [ $(type -t run_or_exit) != "function" ]; then
   echo Could not find function run_or_exit > /dev/stdout
   exit 1
fi

input="${1%.d}.d" 
output="${1%.d}.xvg" 
echo Convert $input to $output
n_lines=$(wc -l $input | awk '{print 5*($1-1)}')
echo Spline lines are $n_lines for ${atom1}-${atom2}
spline -n $n_lines $input > smooth_${input}
run_or_exit --log log_table_to_xvg_${atom1}_${atom2} $table_to_xvg smooth_${input} $output || exit 1

