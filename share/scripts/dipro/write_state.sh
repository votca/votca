#!/bin/bash
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
cd frame$2

for pair in `cat pairlist`;
do
# parse XML-file
    
    #holes
    if [ "$1" = "h" ]; then
	J=`cat transfer_integrals/$pair.xml | grep "<J>" | head -n 1 | sed -e 's/<J>//g' | sed -e 's,</J>,,g'`
    elif [ "$1" = "hdeg" ]; then
	Jsq=`cat transfer_integrals/$pair.xml | grep "<J_sq_boltz>" | head -n 1 | sed -e 's/<J_sq_boltz>//g' | sed -e 's,</J_sq_boltz>,,g'`
	J=`echo $Jsq | awk '{print sqrt($1)}'`
    elif [  "$1" = "e" ]; then
    #electrons
	J=`cat transfer_integrals/$pair.xml | grep "<J>" | tail -n 1 | sed -e 's/<J>//g' | sed -e 's,</J>,,g'`
    elif [ "$1" = "edeg" ]; then
	Jsq=`cat transfer_integrals/$pair.xml | grep "<J_sq_boltz>" | tail -n 1 | sed -e 's/<J_sq_boltz>//g' | sed -e 's,</J_sq_boltz>,,g'`
	J=`echo $Jsq | awk '{print sqrt($1)}'`
    else
	echo "Invalid option for charge carrier TYPE specified. Aborting!"
    fi

    #get segment IDs
    tmp=`echo $pair | sed -e 's/pair_//g'`
    seg1=${tmp%_*}
    seg2=${tmp#*_}

    #write to statefile

    #get pair_id
    pair_ID=`sqlite3 ../$3 "SELECT _id FROM pairs WHERE conjseg1=$seg1 AND conjseg2=$seg2"`

    # if pair_ID empty, STOP
    if [ ! "$pair_ID" = "" ]; then

       #DELETE OLD ENTRY
       sqlite3 ../$3 "DELETE FROM pair_integrals WHERE pair=$pair_ID"

       #INSERT transfer integral
       sqlite3 ../$3 "INSERT INTO pair_integrals (pair,num,J) VALUES ($pair_ID,0,$J)"
    else
	echo "Empty pair_ID for $pair? ABORTING"
	exit
    fi


done

cd ..