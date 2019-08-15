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

cp $VOTCASHARE/xtp/xml/mapchecker.xml OPTIONFILES/
xtp_run -e mapchecker -o OPTIONFILES/mapchecker.xml -f state.hdf5


cp $VOTCASHARE/xtp/xml/eanalyze.xml OPTIONFILES/
changeoption states "e h" OPTIONFILES/eanalyze.xml
changeoption resolution_sites 0.03 OPTIONFILES/eanalyze.xml
xtp_run -e eanalyze -o OPTIONFILES/eanalyze.xml -f state.hdf5

#running ianalyze

cp $VOTCASHARE/xtp/xml/ianalyze.xml OPTIONFILES/
changeoption states "e h" OPTIONFILES/ianalyze.xml
changeoption resolution_logJ2 0.1 OPTIONFILES/ianalyze.xml
xtp_run -e ianalyze -o OPTIONFILES/ianalyze.xml -f state.hdf5

#copy kmcmultiple from votca share folder to here

cp $VOTCASHARE/xtp/xml/kmcmultiple.xml OPTIONFILES/

changeoption runtime  1000 OPTIONFILES/kmcmultiple.xml
changeoption outputtime  10 OPTIONFILES/kmcmultiple.xml
changeoption field  "10 0 0" OPTIONFILES/kmcmultiple.xml
changeoption carriertype  "electron" OPTIONFILES/kmcmultiple.xml
#run kmcmultiple calculator

xtp_run -e kmcmultiple -o OPTIONFILES/kmcmultiple.xml -f state.hdf5




