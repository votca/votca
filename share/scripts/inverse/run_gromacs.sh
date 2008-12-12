#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script runs gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo Needs: run_or_exit, mdrun
   exit 0
fi

for exe in mdrun; do
   if [ -z $(type -p $exe) ]; then
      echo Could not find $exe > /dev/stderr
      exit 1
   fi
done

if [ $(type -t run_or_exit) != "function" ]; then
   echo Could not find function run_or_exit > /dev/stderr
   exit 1
fi

run_or_exit mdrun
