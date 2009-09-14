#!/bin/bash

if [ "$1" = "--help" ]; then
  echo Functions useful for gromacs 4.0
  echo NEEDS:
  #we add \$GMXDATA here, because gromacs will need it
  echo USES: sed die \$GMXDATA
  echo PROVIDES: get_from_mdp
  exit 0
fi 

check_deps $0

get_from_mdp() {
  local res
  [[ -n "$1" ]] || { echo What?; exit 1;}
  res=$(sed -n -e "s#[[:space:]]*$1[[:space:]]*=[[:space:]]*\(.*\)\$#\1#p" grompp.mdp | sed -e 's#;.*##') || die "get_from_mdp failed" 
  [[ -n "$res" ]] || die "get_from_mdp: could not fetch $1"
  echo "$res"
}

export -f get_from_mdp
