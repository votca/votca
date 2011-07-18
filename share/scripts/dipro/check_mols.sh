#!/bin/bash

#for j in {0..$1};
for (( j=$1; j<=$2; j++ ))
do
    echo "Checking mol_$j"
    if [ "$3" == "T" ]; then
	if [ -e mol_$j/mos ];then
	    test=$( head mol_$j/mos | grep 'scfconv' )
	    if [ "$?" -eq "1" ]; then
		echo "mol_$j" >> TROUBLE.mol
	    fi
	    cat mol_$j/mos > /dev/null 
	    if [ "$?" -eq "1" ]; then
		echo "mol_$j" >> TROUBLE.mol
		cd mol_$j
		rm mos control basis auxbasis energy ridft.out statistics
		cd ..
	    fi
	else
	    echo "mol_$j" >> TROUBLE.mol
	fi
    elif [ "$3" == "G" ]; then
	if [ ! -e mol_$j/fort.7 ];then
	    echo "mol_$j" >> TROUBLE.mol
	fi
    fi
done

