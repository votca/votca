#!/bin/bash
rm -r energies
mkdir energies
cd energies

#set equi to 0
equi=0

#see if $1 is set, set it to first argument
if [[ -n $1 ]]; then
equi=$1
fi


echo "Bond" | gmx energy -f ../ener.edr -b $equi -o Bond.xvg
echo "Angle" | gmx energy -f ../ener.edr -b $equi -o Angle.xvg
echo "Ryckaert-Bell." | gmx energy -f ../ener.edr -b $equi -o Ryckaert-Bell.xvg
echo "LJ-14" | gmx energy -f ../ener.edr -b $equi -o LJ-14.xvg
echo "Coulomb-14" | gmx energy -f ../ener.edr -b $equi -o Coulomb-14.xvg
echo "LJ-(SR)" | gmx energy -f ../ener.edr -b $equi -o LJ-\(SR\).xvg
echo "Disper.-corr." | gmx energy -f ../ener.edr -b $equi -o Disper.-corr.xvg
echo "Coulomb-(SR)" | gmx energy -f ../ener.edr -b $equi -o Coulomb-\(SR\).xvg
echo "Coul.-recip." | gmx energy -f ../ener.edr -b $equi -o Coul.-recip.xvg
echo "Potential" | gmx energy -f ../ener.edr -b $equi -o Potential.xvg
echo "Kinetic-En." | gmx energy -f ../ener.edr -b $equi -o Kinetic-En..xvg
echo "Total-Energy" | gmx energy -f ../ener.edr -b $equi -o Total-Energy.xvg
echo "Conserved-En." | gmx energy -f ../ener.edr -b $equi -o Conserved-En.xvg
echo "Temperature" | gmx energy -f ../ener.edr -b $equi -o Temperature.xvg
echo "Pres.-DC" | gmx energy -f ../ener.edr -b $equi -o Pres.-DC.xvg
echo "Pressure" | gmx energy -f ../ener.edr -b $equi -o Pressure.xvg


cd ..
