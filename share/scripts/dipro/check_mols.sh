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

