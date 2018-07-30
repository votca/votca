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
changeoption mode optimize dftgwbse.xml
changeoption molecule CO.xyz dftgwbse.xml
changeoption openmp 4 mbgft.xml
changeoption threads 4 xtpdft.xml
changeoption shift_type iterate mbgft.xml
echo "Running dft + gwbse, output can be found in dftgwbse.log"
xtp_tools -e dftgwbse -o dftgwbse.xml > dftgwbse.log


