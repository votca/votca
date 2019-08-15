#!/bin/bash

#convienience function to change xml option
changeoption(){
    sed -i "s&<${1}.*>.*</${1}>&<${1}>${2}</${1}>&" $3
}
#convienience function to delete xml option
deleteoption(){
 sed -i "s&<${1}.*>.*</${1}>&&" $2
}

echo $VOTCASHARE

#make OPTIONFILE folder, you can put all options into a single options.xml file but experience has shown, that this is better.

mkdir OPTIONFILES

#copy kmcmultiple from votca share folder to here

cp $VOTCASHARE/xtp/xml/kmcmultiple.xml OPTIONFILES/

changeoption runtime  1000 OPTIONFILES/kmcmultiple.xml
changeoption outputtime  10 OPTIONFILES/kmcmultiple.xml
changeoption field  "0 0 0" OPTIONFILES/kmcmultiple.xml
changeoption carriertype  "singlet" OPTIONFILES/kmcmultiple.xml
changeoption rates  "calculate" OPTIONFILES/kmcmultiple.xml
#run kmcmultiple calculator

xtp_run -e kmcmultiple -o OPTIONFILES/kmcmultiple.xml -f state.hdf5

#copy kmclifetime from votca share folder to here

cp $VOTCASHARE/xtp/xml/kmclifetime.xml OPTIONFILES/

changeoption numberofinsertions 5 OPTIONFILES/kmclifetime.xml

#run kmclifetime calculator
# the lifetimes.xml file has to be created from hand.

xtp_run -e kmclifetime -o OPTIONFILES/kmclifetime.xml -f state.hdf5


