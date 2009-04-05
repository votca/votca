#!/bin/bash

get_from_mdp() {
   [[ -n "$1" ]] || { echo What?; exit 1;}
   sed -n -e "s#[[:space:]]*$1[[:space:]]*=[[:space:]]*\(.*\)\$#\1#p" grompp.mdp | sed -e 's#;.*##' || exit 1
   exit 0
}

export -f get_from_mdp
