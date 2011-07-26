#!/bin/bash

cd frame$2

for pair in `cat pairlist`;
do
# parse XML-file
    
    #holes
    if [ "$1" = "h" ]; then
	J=`cat transfer_integrals/$pair.xml | grep "<J>" | head -n 1 | sed -e 's/<J>//g' | sed -e 's,</J>,,g'`
    elif [  "$1" = "e" ]; then
    #electrons
	J=`cat transfer_integrals/$pair.xml | grep "<J>" | tail -n 1 | sed -e 's/<J>//g' | sed -e 's,</J>,,g'`
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