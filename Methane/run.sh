#!/bin/bash

echo $VOTCASHARE

xtp_map -t MD_FILES/topol.tpr -c MD_FILES/conf.gro -s system.xml -f state.sql

mkdir OPTIONFILES

cp $VOTCASHARE/xtp/xml/neighborlist.xml OPTIONFILES/

sed -i "s_<constant>*</constant>_<constant>0.3</constant>_" OPTIONFILES/neighborlist.xml
