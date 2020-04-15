#!/bin/bash

#convienience function to change xml option
changeoption(){
    sed -i "s&<${1}.*>.*</${1}>&<${1}>${2}</${1}>&" $3
}

cp $VOTCASHARE/xtp/xml/dftgwbse.xml .

changeoption molecule methane.xyz dftgwbse.xml
echo "Running dft + gwbse, output can be found in dftgwbse.log"
xtp_tools -e dftgwbse -o dftgwbse.xml> dftgwbse.log


cp $VOTCASHARE/xtp/xml/partialcharges.xml .
cp $VOTCASHARE/xtp/packages/esp2multipole.xml .

changeoption input methane.orb partialcharges.xml
changeoption output methane.mps partialcharges.xml

echo "Running CHELPG fit"
xtp_tools -e partialcharges -o partialcharges.xml

cp $VOTCASHARE/xtp/xml/gencube.xml .

changeoption input methane.orb gencube.xml
changeoption output methane.cube gencube.xml
echo "Running cube file generation"
xtp_tools -e gencube -o gencube.xml

cp $VOTCASHARE/xtp/xml/spectrum.xml .
changeoption orbitals methane.orb spectrum.xml
changeoption upper 25 spectrum.xml
changeoption lower 9 spectrum.xml
changeoption points 1000 spectrum.xml
echo "running spectrum calculation"
xtp_tools -e spectrum -o spectrum.xml
