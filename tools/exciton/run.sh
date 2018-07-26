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
cp $VOTCASHARE/xtp/packages/mbgft.xml


