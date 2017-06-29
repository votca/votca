#!/bin/bash -e

if [[ -d energies ]]; then
  echo "Removing energies directory"
  rm -rf energies
fi

mkdir energies
cd energies

equi=0
#see if $1 is set, set equi to first argument, other use default 0
[[ -n $1 ]] && equi=$1 && shift

for x in "LJ-(SR)" "Disper.-corr." "Coulomb-(SR)" "Potential" "Kinetic-En." "Total-Energy" "Temperature" "Pres.-DC" "Pressure" "$@"; do
  echo "${x}" | gmx energy -f ../ener.edr -b "$equi" -o "${x}.xvg"
done

cd ..
