#! /bin/bash

if [ "$1" = "--help" ]; then
  echo This is a wrapper to find the right file to a task
  echo Usage: ${0##*/} todo method
  echo Uses: \$CSGSCRIPTDIR \$CSGSHARE
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
    echo "Method '$1' unknown" > /dev/stderr
    exit 1;;
esac
}

post_case() {
case $1 in
  update)
    scriptname="${scriptname}_update";;
  add)
    scriptname="${scriptname}_add";;
  *)
    echo "Post case '$1' unknown" > /dev/stderr
    exit 1;;
esac
}

simulation_case() {
case $1 in
  gromacs)
    scriptname="${scriptname}_gromacs";;
  *)
    echo "Simlation programm '$1' unknown" > /dev/stderr
    exit 1;;
esac
}

scriptname=""
called_direct="no"
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
  add_pot)
    scriptname="add_pot"
    shift
    method_case $1;;
  run)
    scriptname="run"
    shift
    simulation_case $1;;
  post)
    scriptname="post"
    shift
    post_case $1;;
  postupd|postadd)
    #no check
    if [ -n "$2" ]; then
      scriptname="$1_$2"
    else
      echo Error: Missing scriptname > /dev/stderr
      exit 1
    fi
    shift 2;;
  functions)
    scriptname="functions"
    shift
    simulation_case $1;;
  --direct)
    called_direct="yes"
    if [ -n "$2" ]; then
      scriptname="$2"
    else
      echo Error: Missing scriptname > /dev/stderr
      exit 1
    fi
    shift 2;;
  *)
    echo "Function '$1' unknown" > /dev/stderr
    exit 1;;
esac

#add ending .sh, but not twice 
[[ "$called_direct" = "yes" ]] || scriptname="${scriptname%.sh}.sh"

#first local CSGSCRIPTDIR. then CSGSHARE
if [ -n  "${CSGSCRIPTDIR}" ] && [ -f "${CSGSCRIPTDIR}/${scriptname}" ]; then
  echo ${CSGSCRIPTDIR}/${scriptname}
  exit 0
elif [ -n  "${CSGSHARE}" ] && [ -f "${CSGSHARE}/${scriptname}" ]; then
  echo ${CSGSHARE}/${scriptname}
  exit 0
else
  echo ${0##*/}: Could not find script ${scriptname} > /dev/stderr
  exit 1
fi
