#!/bin/bash

if [ "$1" = "--help" ]; then
  #we add \$GMXDATA in USES, because gromacs will need it
cat <<EOF
${0##*/}, version @version@
Functions useful for gromacs 4.0

NEEDS:

USES: sed die \$GMXDATA

PROVIDES: get_from_mdp
EOF
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
