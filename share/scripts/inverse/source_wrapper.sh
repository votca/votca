#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This is a wrapper to source the right file
   echo Usage: ${0##*/} todo method
   echo Uses: \$scriptdir
   exit 0
fi

method_case() {
case $1 in
   ibm)
      scriptname="${scriptname}_ibm";;
   imc)
      scriptname="${scriptname}_imc";;
   own)
      scriptname="${scriptname}_own";;
   gromacs)
      scriptname="${scriptname}_gromacs";;
   *)
      echo Method \"$1\" unknown > /dev/stderr
      exit 1;;
esac
}

simulation_case() {
case $1 in
   gromacs)
      scriptname="${scriptname}_gromacs";;
   *)
      echo Simlation programm \"$1\" unknown > /dev/stderr
      exit 1;;
esac
}

scriptname=""
case $1 in
   init*)
      scriptname="initalize"
      shift
      method_case $1;;
   update)
      scriptname="update"
      shift
      method_case $1;;
   prepare)
      scriptname="prepare"
      shift
      simulation_case $1;;
   rdf)
      scriptname="calc_rdf"
      shift
      simulation_case $1;;
   pressure)
      scriptname="calc_pressure"
      shift
      simulation_case $1;;
   convert_potential)
      scriptname="potential_to"
      shift
      simulation_case $1;;
   run)
      scriptname="run"
      shift
      simulation_case $1;;
   *)
      echo Function \"$1\" unknown > /dev/stderr
      exit 1;;
esac


if [ -f "${scriptdir}/${scriptname}.sh" ]; then
   echo ${scriptdir}/${scriptname}.sh
   exit 0
elif [ -f "${scriptname}.sh" ]; then
   echo ${PWD}/${scriptname}.sh
   exit 0
elif [ -f "${CSGSHARE}/${scriptname}.sh" ]; then
   echo ${CSGSHARE}/${scriptname}.sh
   exit 0
else
   echo Could not find script ${scriptname}.sh > /dev/stderr
   exit 1
fi