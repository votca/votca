#!/bin/bash

#convienience function to change xml option
changeoption(){
    sed -i "s&<${1}.*>.*</${1}>&<${1}>${2}</${1}>&" $3
}
#convienience function to delete xml option
deleteoption(){
 sed -i "s&<${1}.*>.*</${1}>&&" $2
}


cp $VOTCASHARE/xtp/xml/dftgwbse.xml .
cp $VOTCASHARE/xtp/packages/xtpdft.xml .
cp $VOTCASHARE/xtp/packages/mbgft.xml .

changeoption dftpackage xtpdft.xml dftgwbse.xml
changeoption molecule methane.xyz dftgwbse.xml
changeoption dftlog system_dft.orb dftgwbse.xml
changeoption openmp 1 mbgft.xml
changeoption threads 1 xtpdft.xml
changeoption shift_type fixed mbgft.xml
echo "Running dft + gwbse, output can be found in dftgwbse.log"
xtp_tools -e dftgwbse -o dftgwbse.xml > dftgwbse.log


cp $VOTCASHARE/xtp/xml/partialcharges.xml .
cp $VOTCASHARE/xtp/packages/esp2multipole.xml .

changeoption openmp 1 esp2multipole.xml
changeoption output methane.mps partialcharges.xml

echo "Running CHELPG fit" 
xtp_tools -e partialcharges -o partialcharges.xml 

cp $VOTCASHARE/xtp/xml/gencube.xml .

changeoption output methane.cube gencube.xml

xtp_tools -e gencube -o gencube.xml

