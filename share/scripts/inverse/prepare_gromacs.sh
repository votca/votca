#! /bin/bash

if [ "$1" = "--help" ]; then
  echo This script implemtents the function initialize
  echo for the Inverse Boltzmann Method
  echo Usage: ${0##*/} $last_sim_dir
  echo Needs: \$atoms, run_or_exit, make_ndx, grompp,sort,uniq
  exit 0
fi

if [ -z "$1" ]; then
  echo Missing argument for ${0##*/} > /dev/stderr
  exit 1
fi

for exe in make_ndx grompp; do
   if [ -z $(type -p $exe) ]; then
      echo Could not find $exe > /dev/stderr
      exit 1
   fi
done

if [ $(type -t run_or_exit) != "function" ]; then
   echo Could not find function run_or_exit > /dev/stderr
   exit 1
fi

cp ../${1}/confout.gro ./conf.gro || exit 1

#realy hacky but it works, looking for a smarter way
atoms="$(for_all non-bonded "echo -e \"\${type1}\\n\${type2}\"" | sort | uniq)" || exit 1
echo -e "a ${atoms}\nq" | make_ndx -f conf.gro &> log_make_ndx || exit 1

run_or_exit grompp -v -n index.ndx || exit 1
