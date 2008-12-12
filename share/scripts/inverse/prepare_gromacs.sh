#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function initialize
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo Needs: \$atoms, run_or_exit, make_ndx, grompp, \$last_dir
   exit 0
fi

for exe in make_ndx grompp; do
   if [ -z $(type -p $exe) ]; then
      echo Could not find $exe > /dev/stderr
      exit 1
   fi
done

if [ -z "$scriptdir" ]; then
   echo scriptdir not defined > /dev/stderr
   exit 1
fi

if [ -z "$last_dir" ]; then
   echo last_dir not defined > /dev/stderr
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

cp ../$last_dir/confout.gro ./conf.gro
echo -e "a ${atoms[*]}\nq" | make_ndx -f conf.gro &> log_make_ndx

run_or_exit grompp -v -n index.ndx
