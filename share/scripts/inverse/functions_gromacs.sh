#!/bin/bash

get_from_mdp() {
  local res
  [[ -n "$1" ]] || { echo What?; exit 1;}
  res=$(sed -n -e "s#[[:space:]]*$1[[:space:]]*=[[:space:]]*\(.*\)\$#\1#p" grompp.mdp | sed -e 's#;.*##') || die "get_from_mdp failed" 
  [[ -n "$res" ]] || die "get_from_mdp: could not fetch $1"
  echo "$res"
}

check_for $0 \$GMXDATA
export -f get_from_mdp
