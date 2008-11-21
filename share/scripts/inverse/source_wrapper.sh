#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This is a wrapper to source the right file
   echo Usage: ${0##*/} todo method
   exit 0
fi

method_case() {
case $1 in
   ibm)
      to_source="${to_source}_ibm";;
   imc)
      to_source="${to_source}_imc";;
   own)
      to_source="${to_source}_own";;
   *)
      echo Method \"$1\" unknown > /dev/stderr
      exit 1;;
esac
}

simulation_case() {
case $1 in
   gromacs)
      to_source="${to_source}_gromacs";;
   *)
      echo Simlation programm \"$1\" unknown > /dev/stderr
      exit 1;;
esac
}

to_source=""
case $1 in
   init*)
      to_source="initalize"
      shift
      method_case $1;;
   update)
      to_source="update"
      shift
      method_case $1;;
   convert_potential)
      to_source="potential_to"
      shift
      simulation_case $1;;
   run)
      to_source="run"
      shift
      simulation_case $1;;
   *)
      echo Function \"$1\" unknown > /dev/stderr
      exit 1;;
esac


if [ -f "${to_source}.sh" ]; then
   echo ${to_source}.sh
   exit 0
elif [ -f "${CSGSHARE}/${to_source}.sh" ]; then
   echo ${CSGSHARE}/${to_source}.sh
   exit 0
else
   echo Could not find script ${to_source}.sh > /dev/stderr
   exit 1
fi