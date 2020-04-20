#!/bin/bash

#convienience function to change xml option
changeoption(){
    sed -i "s&<${1}.*>.*</${1}>&<${1}>${2}</${1}>&" $3
}

echo "Running dft + gwbse, output can be found in dftgwbse.log"
cp $VOTCASHARE/xtp/xml/dftgwbse.xml .
xtp_tools -n methane -e dftgwbse -o dftgwbse.xml> dftgwbse.log

echo "Running CHELPG fit"
cp $VOTCASHARE/xtp/xml/partialcharges.xml .
cp $VOTCASHARE/xtp/packages/esp2multipole.xml .
xtp_tools -n methane -e partialcharges -o partialcharges.xml

echo "Running cube file generation"
cp $VOTCASHARE/xtp/xml/gencube.xml .
xtp_tools -n methane -e gencube -o gencube.xml

echo "running spectrum calculation"
cp $VOTCASHARE/xtp/xml/spectrum.xml .
changeoption upper 25 spectrum.xml
changeoption lower 9 spectrum.xml
changeoption points 1000 spectrum.xml
xtp_tools -n methane -e spectrum -o spectrum.xml
