#!/bin/bash

#convienience function to change xml option
changeoption(){
    sed -i "s&<${1}.*>.*</${1}>&<${1}>${2}</${1}>&" $3
}

echo "Running dft + gwbse, output can be found in dftgwbse.log"
cp $VOTCASHARE/xtp/xml/dftgwbse.xml .

changeoption ranges full dftgwbse.xml
# Change enery mode to optimize
xtp_tools -n CO -e dftgwbse -o dftgwbse.xml -t 4 > dftgwbse.log


